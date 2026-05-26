# 计算结构生物学/分子动力学实验最佳实践

> 基于 `naked-mole-rat-cgas-trim41-simulation` (2026-04-22~05-08) 和 `cjc-1295` (2026-05-13~17) 两个项目的经验总结。
> 累计模拟规模约 5.5 μs，涵盖蛋白-蛋白对接、MD 模拟、FEP、伞形采样、Rosetta 折叠等。

---

## 一、工具选择与性能

### 1. OpenMM vs GROMACS：选引擎

- **OpenMM 比 GROMACS 快 30-50%。** 在相同体系（~324k atoms, Amber14SB + TIP3P, RTX 3090）上实测：OpenMM 44 ns/day vs GROMACS 33 ns/day。原因包括单 GPU 设计、mixed precision、无 domain decomposition 开销。
- **OpenMM 的 Python API 更灵活。** 自定义 restraint（flat-bottom COM, harmonic distance）、CustomCentroidBondForce 用于 umbrella sampling 等场景，OpenMM 表达力远超 GROMACS mdp 文件。
- **GROMACS 生态更成熟。** `gmx` 命令行工具体系完备（rms, distance, hbond, trjconv, bar 等），分析工具链比 OpenMM + MDAnalysis 组合更一体化。FEP 的 dual-topology 在 GROMACS 中直接可用。
- **推荐**：生产 MD 优先 OpenMM（快）；需要 FEP、复杂分析时用 GROMACS；跨引擎验证两者都跑。

### 2. 多 GPU 并行不要 pin CPU

- **`gmx mdrun -pin on` 在多 GPU 并行场景下是性能陷阱。** CJC-1295 项目实测：4 个 mdrun 进程同时跑在 4 张 RTX 3090 上，`-pin on -ntomp 16` 导致 GPU 利用率仅 53%，性能 19 ns/day。去掉 `-pin on` + 降低 `-ntomp` 到 8 后，GPU 利用率恢复到 93%，性能 33 ns/day (+73%)。
- 根因：固定线程到特定 CPU 核心导致多进程间资源竞争。
- **GPU 利用率是最直接的诊断指标**：<60% 几乎肯定意味着 CPU 瓶颈或竞争。

### 3. GROMACS 性能调优清单

- `nstlist` 越大越快，但需配合 `rlist`（如 `nstlist=40, rlist=1.1` 或 `nstlist=50, rlist=1.2`）
- `ntmpi=1` 在单 GPU 场景下通常最优
- 不要在 checkpoint 续跑时尝试改 `dt`——从 checkpoint 提取帧后新建 tpr 会丢失速度信息，导致性能反而下降。正确做法是 `gmx convert-tpr`。
- NHMR (hydrogen mass repartitioning) 允许 dt=0.004，但需提前做，不能从已有 checkpoint 改。

### 4. 硬件选择

- **RTX 3090 (24GB) 性价比极高。** 116k atoms 体系 ~152 ns/day（OpenMM），~320k atoms 体系 ~33 ns/day（GROMACS）。
- Apple Silicon (M3 Pro) 的 OpenCL 后端比 CUDA 慢 3-4×，且 `mixed` precision 不可用，仅适合小规模测试。
- 116k atoms 体系在 RTX 3090 上仅占 ~2-3 GB 显存，24 GB 绰绰有余。

---

## 二、力场与非标准残基

### 5. 非标准氨基酸力场处理

- **D-氨基酸**：
  - **GROMACS + Amber14SB**：标准 `amber14sb.ff` 不含 D-氨基酸模板。需要手动在 `.rtp` 中添加 `[ DALA ]` 条目（复制 `[ ALA ]`），并在 `aminoacids.hdb` 中添加氢原子规则。L/D 的区别仅在于初始坐标——Amber 的 improper dihedral 会维持初始手性，无需单独力场参数。
  - **注意**：GROMACS 2026 conda 包的 `amber14sb.ff/aminoacids.hdb` 存在上游 bug（ALA 条目重复、CSER 混入 THR 规则），需用 `amber19sb.ff` 的干净版本替换。
  - **OpenMM**：标准 `amber/protein.ff14SB.xml` 不含 D-氨基酸模板；`charmm/charmm36_protein_d.xml` 有 DALA 但需改残基名。
  - **Rosetta**：原生支持，将 PDB 残基名 `ALA` 改为 `DAL`，自动应用 `D_AA` chiral_flip patch。
  - **AmberTools (tleap)**：标准 `leaprc.protein.ff14SB` 不含 D-ALA。需使用 `Glycam` 或手动添加。
  - **验证方法**：用 signed_volume `(N−CA)·[(C−CA)×(CB−CA)]` 作为 Cα 手性独立判据——正值=L-Ala，负值=D-Ala。

- **磷酸化残基（SEP）**：
  - AMBER 的 `leaprc.phosaa19SB` 原生支持磷酸化氨基酸。
  - 磷酸化引入 -2 电荷，可能导致静电过度估计。用 S→E 突变（-1 电荷）作为电荷效应对照。
  - MM-GBSA 计算时需使用 protein-only 拓扑（strip 水/离子），因为 sander 不支持 OPC 水的 EP 原子类型。

- **异肽键 (isopeptide bond)**：
  - E2~Ub 复合物中 Ub C-ter 与 E2 活性位点 Cys 的异肽键在标准 tleap/pdb2gmx 中不会被识别。
  - 需要用 harmonic distance restraint 强制维持 ~1.6 Å 的共价键距离（在 OpenMM 中用 CustomBondForce）。

### 6. PDB 格式是永恒的痛点

- **Bio.PDB 的列对齐问题**：Python Bio.PDB 的 PDB writer 在多项目中出现至少 3 次列对齐错误（残基名变成 YS 而非 LYS、chain ID 丢失、residue ID 与原始不同）。每次 debug 耗时 30-60 分钟。
- **建议**：项目早期写一个经过列验证的 PDB writer/validator，后续所有脚本统一调用。一次性投入 1 小时，省去反复 debug 的时间。
- **pdb4amber 后务必验证** chain ID 和 residue ID 是否正确传递。
- **Parmed 转换 GROMACS→OpenMM 时的 PDB wrap 问题**：`gmx editconf` 输出 PDB 时坐标被包装（wrapped），导致 OpenMM 键能爆炸（34.8 亿 kJ/mol）。必须用 `gmx trjconv -pbc whole` 先恢复分子完整性。

---

## 三、多引擎交叉验证

### 7. 跨引擎验证是项目救星

- **GROMACS vs OpenMM 的不一致**引导发现了两个严重 bug：
  1. **CMAP bug**：parmed 将 ff19SB 的 14 种残基特异性 CMAP 压缩为 1 种，GROMACS RMSD 被高估 4×
  2. **Python 对齐 bug**：`rotation_matrix()` 在未中心化坐标上计算，RMSD 被高估 82%
- 没有跨引擎验证这一步，大量数据会被错误使用。
- **推荐做法**：用两种不同引擎跑同一个体系 10-20ns，确认 RMSD/COM/Rg 在 20% 以内一致，再投入大规模生产。

### 8. 跨引擎验证中的常见陷阱

- **力场转换工具可能引入系统性错误。** parmed 的 AMBER→GROMACS 转换丢失了残基特异性 CMAP。始终用原生实现（GROMACS 2026 内置 `amber19sb.ff` 或 `amber14sb.ff`）作为 gold standard。
- **NPT vs NVT**：两引擎的 production ensemble 必须一致。OpenMM 用 NVT（移除 barostat），GROMACS 也应用 `pcoupl = no`。
- **约束精度**：GROMACS 的默认 `lincs_iter=1, lincs_order=4` 对长模拟可能不够，建议 `lincs_iter=2, lincs_order=6`。
- **分析脚本必须用独立方法验证。** 在用于正式结论前，Python 写的 RMSD/alignment 要与原生工具（`gmx rms`, MDAnalysis `alignto()`）核对。单一工具的分析结果不可信。

### 9. GROMACS 转化 OpenMM 必须处理 PBC 包裹

- **GROMACS 默认输出 wrapped 坐标，OpenMM DCD 输出 unwrapped 坐标。**
- 未修复 PBC 前，GROMACS COM 被误报为 28 Å（真实 44.6 Å），被误判为"RMSD ratio 5.8× 的严重偏离"。
- **解决**：`gmx trjconv -pbc mol` 或 `MDAnalysis.transformations.unwrap`。
- **警示**：分析 GROMACS 轨迹时，**必须先修复 PBC**，否则所有基于 COM/Rg/RMSD 的结论都是错的。

### 10. OpenMM 复刻验证的最佳实践

- 同力场、不同引擎、独立实现的复刻是 paper-grade cross-validation。温度、密度、能量在同一水平，可有效反驳"XX 引擎数值伪迹"的质疑。
- 仅验证热化学层面（温度/密度）不够——需要把复刻跑到至少 10 ns，看关键几何指标（如催化距离）的时序是否与主引擎重合。
- OpenMM 的 MonteCarloBarostat 与 GROMACS 的 Parrinello-Rahman 会给出 ~2% 的密度差异，这是系统误差，非 bug。

---

## 四、实验设计

### 11. Replica 设计是不可妥协的

- **S305E 的教训**：3 个 replica 分别显示完全解离（ΔG = +7.2）、稳定结合（ΔG = -22.9）、弱结合（ΔG = -11.9）。单 replica 会给出完全误导的结论。
- 每个条件至少 3 replica，并报告 replica 间的变异系数。
- 不要用独立样本 t-test 处理自相关的 MD 轨迹——用 correlated t-test 或有效样本量校正。

### 12. 渐进式模型构建

- 不要一步到位构建最复杂的系统。逐层加复杂度：
  1. 二元复合物 → 验证界面稳定性
  2. + RING + 异肽键 → 验证 E2~Ub 稳定性
  3. + SPRY + CC linker → 验证锚定效应
  4. US/PMF → 定量自由能面
- 每一层回答前一层的问题，同时暴露下一层需要的改进。
- **先做 frame-0 验证再投 GPU 时间**：复杂结构拼装后，先做 clash check + 关键距离验证。一次 Kabsch 刚体对齐的灾难性 clash（999+ contacts）在 frame-0 就能发现，而不是浪费半天 GPU 后才发现。

### 13. Rosetta 折叠/对接的采样量

- **FlexPepDock 的 ab-initio 模式对长肽不可靠。** 29 残基 GHRH 在 960 个 decoy 中没有一个保持催化几何——即使起始结构是正确的。FlexPepDock 的 low-res preoptimize 阶段会大幅扰动初始位置。
- **推荐**：对 >10 残基的肽段，不要信任 ab-initio 对接。用晶体结构起始 + constrained refinement（CoordinateConstraint on backbone CA）。
- **全局对接 10 个 decoy 通常已足够**看到主要结合模式。人类 WT 的 Fnat 0.89（高刚性界面）vs Hgal 的 Fnat 0.39（柔性界面）——decoys 之间的 variability 本身就是信息。

### 14. 结构预测工具的交叉验证

- **AF3、Boltz-2、Chai-1 对低亲和力/瞬态复合物的界面预测是定性不同的**（Jaccard = 0.00）。单一模型无法给出可靠结论。
- 对于 ipTM < 0.33 的预测，必须用至少两种方法交叉验证。
- 这些模型对单体结构的预测高度可靠（cGAS RMSD ~1 Å），但对界面几何的预测不可靠——这正是对接 + MD 策略有价值的前提。
- **AF3 突变序列验证**：4 个点突变对整体 backbone 几何几乎无影响（ΔRMSD < 0.3 Å）——这在实验前不知道，但 AF3 直接排除了"突变改变活性位点几何"这个假说方向。

### 15. 不要信任序列编号

- **PDB / FASTA / 论文中的残基编号系统可能不同步。** 裸鼹鼠项目初期将人源 4mut 标记为 C463S，经全局序列比对才纠正为 D431S（偏移 32 位）。
- 项目第一天就应该做全局序列比对确认编号。错误的编号会导致错误的 fasta、错误的脚本、错误的数据解读——Hgal_rev RMSD 从 13.04 Å（基于错误位点）修正为 0.98 Å（正确位点），结论完全不同。
- 建立"编号字典"：论文编号 ↔ 全长编号 ↔ 构建体编号 ↔ topology resid，全局引用。

---

## 五、分析脚本与代码管理

### 16. 分析脚本验证铁律

- **任何分析脚本在用于正式结论前，必须用独立方法（原生工具或手工计算）验证。**
  - Python 写的 RMSD/alignment → 与 `gmx rms` 或 MDAnalysis `alignto()` 核对
  - 自定义 Kabsch 实现 → 与已知基准值对比
  - 距离计算 → 与 PyMOL 手动测量核对
- cGAS-TRIM41 项目中，Python 对齐 bug 导致 GROMACS RMSD 被高估 82%，差一点让 200ns 的正确数据被当"严重偏离"而废弃。

### 17. 共享库必须立即接入

- **lib/ 创建后立即接入至少 2-3 个实际脚本**，否则永远是"下周我会用"的僵尸代码。
- cGAS-TRIM41 项目中整理了 `lib/`（paths.py, stats.py, mda_tools.py, plot_style.py），但整理后发现：**没有一个生产脚本实际 import 它**。33 个脚本仍然各自硬编码 `/home/scroll/...` 路径。
- 建库的同时写 CI 检查（或至少 grepp 确认没有硬编码路径逃逸）。

### 18. 脚本组织

- 按功能分目录：`01_build/` → `02_md/` → `03_analysis/` → `04_viz/` → `05_utils/`
- 旧版本归档到 `archive/_versions/`，不要留在根目录。64 个 v1-v10 脚本比 41 个干净脚本难维护得多。
- 实验脚本（一次性的 mock/dry-run/test）归档到 `archive/_experiments/`。

### 19. 自动化调度要谨慎

- **auto-launcher 存在多个非平凡陷阱**：
  1. 进程检测逻辑不完整（只检测 `run_production.py` 不检测 `run_md.py`）→ 误判 GPU 空闲 → 重复启动
  2. CUDA context 初始化有 race condition（新进程 ~5-10s 后才出现在 `nvidia-smi pmon` 中）
  3. Python stdout 在 nohup 环境下被块缓冲 → log 文件 1 小时无更新
  4. Launcher 无法处理构建错误（NaN 崩溃后不会自动修复）
- **小规模实验（<10 个 rep）手动管理更可靠。** 自动化调度的 ROI 在大规模（>20 个独立任务）时才为正。

### 20. 文件格式的鲁棒处理

- **DCD 文件头损坏恢复**：前 8 字节为 0x00（应为 CORD 魔数），但帧体完整时可以通过二进制扫描帧结构恢复数据。Hgal_WT rep2 因此挽救了 1008 帧（~100 ns）。
- **gmx trjconv 坐标缩放 bug**：某次 DCD 修复后坐标被放大 10 倍且丢失 box 信息——修复后需用独立的坐标范围 / RMSD / COM 检查确认数据完整性。
- **gmx make_ndx 使用 PDB 残基编号而非序列残基索引**：当使用 `-merge all` 合并多链时，链 A 保留 39-766，链 B 保留 1-29。`r 728` 选中的是 DPP-IV Val728 而非 GHRH Tyr1。**始终用 `gmx dump` 验证 atom indices。**

### 21. GROMACS 原子选择教训

- `gmx make_ndx r N` 选择的是 **PDB residue number**，不是 sequential residue index。
- 多链合并后（`-merge all`），各链保持原始 PDB 编号。Water 从最大残基号+1 开始。
- **验证方法**：`gmx dump -s topol.tpr | grep -A5 "residue (N)"` 确认编号映射。
- **最佳实践**：为分析建独立 index file，所有关键原子（催化残基 OG、底物 C、盐桥原子等）写入 `.ndx` 并注释说明编号来源。避免每次重写 `r N` 选择。

---

## 六、错误检测与自我纠错

### 22. 自我纠错文化

- 每次发现错误 → 承认 → 修复 → 更新文档的循环必须彻底。cGAS-TRIM41 项目最难能可贵的是失败记录和成功记录同样详细——ClusPro 失败、SDOCK2.0 超时、auto-launcher race condition、NaN 崩溃——全部可追溯。
- 被推翻的旧结论**必须在原条目顶部标注 `> ⚠️ SUPERSEDED YYYY-MM-DD: 见XX修正`**。CJC-1295 项目中 5/15 手性"修正"被 5/17 推翻，但 5/15 的条目原样保留，下一个读者会被误导。

### 23. 独立验证的硬卡点

- CJC-1295 项目的手性验证教训：**χ₁ (N-CA-CB-C) 是侧链旋转角，不是手性判据。** L-Ala 和 D-Ala 都可以有正或负的 χ₁。用 signed_volume `(N−CA)·[(C−CA)×(CB−CA)]` 才是 Cα 手性的正确几何判据。
- 当两次独立解读得出相反结论时（5/13 说 -57°=L vs 5/15 说 -55°=D），必须回到物理第一性原理独立验证。
- **任何涉及"标签互换"的修正都必须用物理测量验证**，而不是基于推理链。一条推理链可能逻辑上自洽但物理上全错。

### 24. 常见数据解读陷阱

| 陷阱 | 表现 | 解决 |
|------|------|------|
| 原子索引选错 | "完美 3.15 Å" 实际测的是 Ser630-OG → Asp3-N 而非 Ala2-C | `gmx dump` 验证索引 |
| PBC 包裹误解 | GROMACS COM "28 Å" 实际是 44.6 Å | `gmx trjconv -pbc mol` |
| 对齐 bug | RMSD 高估 82-104% | 与 `gmx rms` gold standard 核对 |
| t-test 对自相关数据 | p 值不可信 | correlated t-test |
| 手性标签误判 | χ₁ 误当手性判据 | signed_volume 独立验证 |
| 数据 key 名错误 | `com` vs `com_dists`，误报"0 frames" | 先 `np.load` 打印所有 keys |

---

## 七、文档与实验记录

### 25. 项目日志规范

- **按 § 编号记录每一次实验、每一个 bug 修复、每一个待决策事项。** cGAS-TRIM41 项目 3600+ 行的 project_log 使得后来者（包括几天后的自己）可以精确理解每个决策的背景。
- 每个条目包含：
  - 日期和时间
  - 实验目的
  - 方法和参数
  - 结果（定量数据）
  - 结论和下一步
- **失败记录同样详细。** 知道"什么不 work"和知道"什么 work"同样有价值。
- 日志末尾标注维护者和最后更新时间。

### 26. 项目文档结构

```
docs/
  00-project/       # 项目概览、日志、回顾、决策记录
  10-reports/       # 单项实验报告
  15-literature/    # 文献调研
  20-protocols/     # 实验方案
  30-diagnostics/   # 问题诊断
  40-reviews/       # 外部评审
  50-infra/         # 硬件/软件/环境
paper/              # 论文手稿
```

### 27. 论文写作时机

- **论文应该在核心数据齐备后写。** 早期草稿会导致：
  - 框架限制思维（"已经有这个 section 了，不舍得删"）
  - 不断修补僵尸框架
  - 引用全是 PLACEHOLDER
  - 突变编号还是错的
- **推荐顺序**：先写 Methods 和 Introduction（数据无关），Results 最后（数据驱动）。

### 28. 外部评审

- **主动请人做 code/doc review**。cGAS-TRIM41 请了 Kimi 和 Gemini review，发现了 t-test 误用、突变编号未清理、docking pose bias、盐浓度注释不符等致命问题。
- Review 发现的问题逐条验证并修复，记录在 `docs/40-reviews/` 下。
- **Review 可以来自 AI**——关键在于被 review 的意愿和回应速度。

---

## 八、项目管理与资源分配

### 29. 知道什么时候停

- **S305 磷酸化是代价最大的教训**：结论在 200 ns 时已明确（完全解离），但跑了 1200 ns+MM-GBSA。对论文的贡献是 1 小节 side note，不应该是 1200 ns 的投入。
- 明确主线后，side projects 用最小可行实验回答，达到结论就停。
- **资源性价比表**定期自审：

| 实验 | GPU 时间 | 论文贡献 | 性价比 |
|------|---------|---------|--------|
| 四系统 MD | ~2.4 μs | 核心 Table | 高 |
| 磷酸化 MD | ~1.2 μs | 1 小节 | **低** |
| 四元 FULL + US | ~0.6 μs | 核心 | 高 |
| K347 探索 | ~0.3 μs | negative result | 中 |

### 30. 数据完整性审计

- 定期做"数据完整性审计"：每类系统有哪些 replica、哪些已完成、哪些有缺失、哪些分析已跑、哪些未跑。
- cGAS-TRIM41 §51.1 的审计表发现 Hgal NEW 构建零生产轨迹——如果不审计，论文写到一半才发现缺数据。
- 数据缺失标记（如 L-Ala 全长缺失 9.4-66.7 ns 的 57 ns 数据）必须在所有引用处显式说明。

### 31. 磁盘空间管理

- DCD 文件头损坏可能源于磁盘瞬时满。320k atom 体系的单条 200ns 轨迹约 8-10 GB（取决于 sampling rate）。
- **定期检查磁盘剩余空间**（特别是 `/tmp` 和输出盘）。
- 已完成轨迹及时压缩或归档，释放工作盘空间。

---

## 九、统计严谨性

### 32. 时间序列自相关

- **MD 轨迹是高度自相关的时间序列，不能用独立样本 t-test。**
- 错误使用：对 200ns 2000-frames 轨迹直接 t-test → p 值严重低估（假阳性）。
- 正确方法：correlated t-test（有效样本量校正，`n_eff = n_frames / (1 + 2∑τ_lag)`）。
- MM-GBSA 的标准差通常为 7-14 kcal/mol（对蛋白-蛋白体系），远大于看起来的 SEM。始终报告 SD 而非 SEM。

### 33. 多重比较校正

- 541 残基的 ΔRMSF 检验：0 个 Bonferroni 显著。19 个未校正显著在 541 重校正后全部不显著。
- 报告时同时给出未校正和校正后的结果，但结论基于校正后。

### 34. FEP 相关

- BAR 内部误差（±0.31 kJ/mol）不是 sampling 误差。需要 time-block convergence 或 reverse FEP (hysteresis check) 独立验证。
- +0.83 ± 0.31 kJ/mol 是**统计显著**（z = 2.68, p ≈ 0.007）但**生物学可忽略**（~0.32 kT, ~1.4× 亲和力差异）。正确表述是 "statistically detectable but biologically negligible"，而不是 "indistinguishable from zero"。
- FEP 在非催化态系统上测的 ΔΔG，物理意义取决于所测的 ensemble 是否代表真实结合态。如果被测系统 0% productive pose，ΔΔG 是 thermodynamic decoupling 实验，不一定反映催化 context。

---

## 十、具体技术清单

### 35. 结构拼装前的检查清单

- [ ] Clash check（< 2 Å 重原子接触数）
- [ ] 关键距离验证（催化距离、盐桥距离、界面接触）
- [ ] Chain ID 一致性（tleap/pdb2gmx 后是否保持）
- [ ] Residue ID 一致性
- [ ] 能量最小化后能量是否收敛（ΔE 是否 < 预设阈值）
- [ ] 初始总能量是否合理（NaN/10²¹ 意味着 PBC/坐标问题）

### 36. MD 启动前检查清单

- [ ] 力场选择合理（蛋白力场 + 水模型 + 离子参数）
- [ ] 溶剂化缓冲层足够（≥10 Å）
- [ ] 体系电荷中和
- [ ] Heating 阶段温度平稳升至目标
- [ ] NPT 阶段密度/体积收敛
- [ ] 前 100 ps 无 NaN、无能量爆炸
- [ ] Checkpoint 写入正常

### 37. 分析前检查清单

- [ ] PBC 已处理（GROMACS 轨迹）
- [ ] 原子索引已验证（`gmx dump`）
- [ ] 分析脚本已与独立方法核对
- [ ] 了解数据的缺失/缺口
- [ ] 关键的 ndx 文件已保存并注释

### 38. GROMACS↔OpenMM 转换清单

- [ ] `gmx trjconv -pbc whole` 恢复 unwrapped 坐标
- [ ] 在 DPP-IV C-ter 和 GHRH N-ter 之间插入 TER record
- [ ] 残基名修复：SOL→HOH, NA→Na+, CL→Cl-
- [ ] 验证初始 PE 合理（不应 > 10⁷ kJ/mol）
- [ ] 力场文件匹配（amber14-all.xml + tip3p.xml）

---

## 十一、最重要的原则

### 39. 多方法交叉验证是不可替代的

- **结构预测**：AF3 + Boltz-2 + Chai-1
- **对接**：LightDock + Rosetta + ClusPro
- **MD 引擎**：OpenMM + GROMACS
- **结合能**：MM-GBSA + Rosetta I_sc + 催化几何分析
- 三个层面互相印证，才是结论可信度的基础。

### 40. 诚实记录一切

- 失败的实验、推翻的假设、废弃的数据、错误的分析——全部如实记录。
- 3 个月后回来看项目的人（可能就是你）需要知道什么做过、什么不 work、为什么当时选择了那条路。
- 被推翻的旧结论不要删除，标注 SUPERSEDED 并指向修正后的结论。
- 这种诚实度在学术计算项目中很少见，但它是项目长期可信度的唯一保障。

---

*最后更新：2026-05-26*
*来源：naked-mole-rat-cgas-trim41-simulation + cjc-1295 两个项目的完整经验回顾*
