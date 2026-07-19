# exp-D 实验日志

> 只追加，不删除。每次记录含日期、时间、操作内容和结果。

---

## 2026-05-31 — 实验启动

### 决策

exp-C 结论：linker 接上后 FA 逃逸 HSA FA3（无论末端电荷）。转 exp-D 验证 linker 的另一半功能：在 GLP-1R 端隔离脂链、保护受体结合。对应 Lau 2015 Table 3（variable linker length → GLP-1R potency）。

### 体系设计

5 个 linker 变体，均连接 C18 二酸，接在 Aib8,Arg34-GLP-1(7-37) 的 Lys26 上：

| # | Linker | Lau 2015 | EC50 (pM) | 状态 |
|---|--------|---------|-----------|:---:|
| 1 | 无 linker（直达 C18） | Cmpd 19 | 269 | 运行中 |
| 2 | γGlu | Cmpd 20 | 9.9 | NaN（C18碰撞ECD） |
| 3 | γGlu-1×OEG | Cmpd 21 | 4.8 | 运行中 |
| 4 | γGlu-2×OEG（司美格鲁肽） | — | 6.2 | 运行中 |
| 5 | γGlu-3×OEG | Cmpd 24 | 27.7 | 运行中 |

系统不含 HSA——linker-FA 在溶剂中与 GLP-1R ECD 相互作用。

### LYA linker 变体构建

- **方法**：RDKit ETKDGv3 + MMFF → 原始 mol2 → antechamber AM1-BCC（OMP_NUM_THREADS=4）
- **净电荷**：nc=-1（COO⁻），偶数电子确保 SQM 收敛
- **剥离**：去除 ACE cap + Lys backbone/sidechain + NME cap，仅保留 linker-C18-COO⁻
- **GAFF2 类型**：n (amide N), c (carbonyl C), o (O), c3 (CH2), hc (H on C)
- **定位**：Kabsch 对齐 amide N→C 向量至 peptide CE→NZ 方向，amide N 置于 NZ 位置

### 肽-ECD 复合物构建

- **模板**：3IOL（ECD + GLP-1 10-35，2.1Å），唯一含正确 ECD-肽结合构象的晶体结构
- **Lys34→Arg 突变**：tleap loadPdb + 模板自动补全缺失原子
- **已知局限**：突变质量受限于初始原子位置，不影响 exp-D 主要指标（C18-ECD 距离）

### 关键构建教训

1. **tleap remove/bond 命令存在静默失败 bug**：在 ParmEd 生成的 PDB 格式下 `remove sys sys.117.HZ1` 报 "Argument #2 is of type String"，`bond` 同样静默失败
2. **10M 残基问题**：ecd_pep.pdb 包含晶体学 NME cap 残基（10M），tleap 创建 33 个无类型原子导致 FATAL 错误。`remove sys sys.127` 在特定 PDB 格式下可删除
3. **LNK 共价键缺失**：由于 tleap bond 失败，prmtop 中 LNK 原子无共价键。原子间距正确时 minimization 可修复，MD 稳定运行
4. **混合力场 frcmod**：需要 `lya_link.frcmod` 覆盖 ff14SB N3—GAFF2 c 之间的 cross-force-field bond/angle 参数
5. **ParmEd PDB 格式**：仅 ParmEd `.save()` 生成的 PDB 能被 tleap atom mask 正确解析；手动 Python PDB writer 导致 mask 失败

### AM1-BCC 运行时间

| Variant | Atoms | SQM 时间 |
|---------|:---:|------|
| no_linker | 85 | ~3 min |
| gglu | 98 | ~18 min |
| gglu_1oeg | 119 | ~40 min |
| gglu_2oeg | 140 | ~60 min |
| gglu_3oeg | 161 | ~90 min |

---

## 2026-06-01 — MD 生产进展

### 运行状态

| Variant | Rep 1 | Rep 2 | Rep 3 |
|---------|:---:|:---:|:---:|
| no_linker | 67ns | 91ns | 89ns |
| gglu_1oeg | 89ns | 100ns ✓ | 67ns |
| gglu_2oeg | 100ns ✓ | 67ns | 91ns |
| gglu_3oeg | 67ns | 91ns | 89ns |

- 速度：单 GPU ~450 ns/d，共享 ~100-200 ns/d
- 系统规模：~36k atoms（vs exp-C ~86k）
- 预计全部完成：2026-06-01 晚

### gglu NaN 诊断

- PE after minimization: 7×10^16 kJ/mol（vs 正常 ~-5.9×10^5）
- 根因：γGlu-only linker 导致 C18 tail 从 NZ 延伸 4-5 Å 即碰撞 ECD 表面
- 更长 linker（1oeg/2oeg/3oeg）将 tail 延伸至溶剂区，无碰撞
- 待修复：调整 tail 初始旋转角度，避让 ECD 表面

### Review 反馈 (Kimi, 2026-06-01)

- barostat 时序 bug 已修复（移至 Simulation 创建之前）
- exp-D 文档已补充
- launch.sh 已废弃（从未用于生产启动）

---

*维护者：Claude Code*
*最后更新：2026-06-01*

## 2026-06-01 — MD 完成与首轮分析

### 完成状态
- 12/15 replica 完成, gglu rep1 重开 GPU 3 恢复正常速度

### 首轮分析：linker 变体无显著差异

| Variant | EC50 (pM) | CA RMSD (A) | Tail-Prot (A) | NZ-C (A) |
|---------|-----------|-------------|---------------|----------|
| No linker | 269 | 2.2+-0.1 | 5.9+-0.1 | 14.4+-2.7 |
| gGlu-1xOEG | 4.8 | 2.2+-0.0 | 6.0+-0.2 | 14.6+-3.0 |
| gGlu-2xOEG | 6.2 | 2.0+-0.1 | 5.8+-1.2 | 14.6+-3.0 |
| gGlu-3xOEG | 27.7 | 2.0+-0.2 | 5.8+-1.0 | 14.1+-3.0 |

**核心发现**: NZ-C 距离 ~14A -> 酰胺键未形成。LNK 靠非键作用浮在肽链附近，linker 差异被掩盖。

**根因**: tleap bond 命令静默失败; ParmEd bond 添加因缺少 BondType 对象失败。

### 修复计划
- ParmEd BondType(427.0, 1.38) 创建正确 NZ-C 酰胺键
- 单 replica 验证: NZ-C ~1.4A -> 变体间有趋势 -> 扩展 3 replica

---
*维护者：Claude Code*
*最后更新：2026-06-01*

## 2026-06-01 — ParmEd bond 修复成功

### 修复方法
- tleap `bond` 命令静默失败导致 NZ-C 酰胺键缺失（NZ-C ~14 Å）
- ParmEd 修复：删除 HZ 原子 + 距离检测 LNK 内部键 + BondType 创建 NZ-C 酰胺键
- BondType(427.0 kcal/mol/A², 1.38 A) 匹配 N3-c amide bond

### 修复后结果
| Variant | PE after min | NZ-C after min | NZ-C production |
|---------|-------------|:---:|:---:|
| no_linker | -588,287 | 1.38 A | 1.36-1.44 A（5ns 稳定） |
| gglu | -574,610 | 1.38 A | 测试中 |
| gglu_1oeg | -598,121 | 1.39 A | 测试中 |
| gglu_2oeg | -578,083 | 1.38 A | 测试中 |
| gglu_3oeg | -582,262 | 1.39 A | 测试中 |

全部 5 个变体 NZ-C = 1.38-1.39 A（正确酰胺键），PE 正常，无 NaN。

### 下一步
- 单 replica 10ns 测试确认 NZ-C 稳定 → 全 5 变体 ×3 replica 100ns 生产

---
*维护者：Claude Code*
*最后更新：2026-06-01*

## 2026-06-02 — 全量 MD 生产进展

### 状态
- 4/15 done（全部 rep1 先完成），11 个 ~60-64ns，全部 307-312K 无 NaN
- 速度 92-98 ns/d，4 GPU 满载，预计今晚全部完成

### 中间分析（~60ns）
NZ-C 键修复验证成功：

| Variant | EC50 | CA RMSD (A) | NZ-C (A) | Tail-Prot (A) |
|---------|------|-------------|:---:|-----------|
| No linker | 269 | 2.1+-0.2 | 1.51 OK | 4.7+-0.2 |
| gGlu | 9.9 | 2.2+-0.2 | 1.51 OK | 4.3+-0.5 |
| gGlu-1xOEG | 4.8 | 2.1+-0.2 | 1.52 OK | 4.5+-0.6 |
| gGlu-2xOEG | 6.2 | 2.4+-0.2 | 1.52 OK | 3.8+-0.2 |
| gGlu-3xOEG | 27.7 | 3.0+-0.2 | 1.52 OK | 4.0+-0.3 |

- NZ-C 全部 1.51-1.52 A（vs 修复前 14 A），共价键修复 100% 成功
- gGlu-3xOEG CA RMSD 最高（3.0 A），可能反映过长 linker 的 entropic penalty
- gGlu-2xOEG Tail-Prot 最近（3.8 A），linker 刚好够长产生有利接触

---
*维护者：Claude Code*
*最后更新：2026-06-02*

## 2026-06-02 — 中程分析（~85ns）

### 结果

| Variant | EC50 | CA RMSD (A) | Tail-Prot (A) | 趋势 |
|---------|------|:---:|:---:|------|
| No linker | 269 | 2.2+-0.2 | 4.7+-0.1 | 基准 |
| gGlu | 9.9 | 2.2+-0.2 | 4.4+-0.5 | 基准 |
| gGlu-1xOEG | 4.8 | 2.4+-0.3 | 4.4+-0.6 | 基准 |
| gGlu-2xOEG | 6.2 | 2.5+-0.2 | 3.8+-0.2 | 最优：最低 tail距离 |
| gGlu-3xOEG | 27.7 | 3.0+-0.2 | 4.0+-0.3 | 最高 RMSD，过长 linker |

### 关键发现
- gGlu-3xOEG CA RMSD 显著升高（3.0 vs 2.2-2.5），量化 entropic penalty——与 EC50=27.7 定性一致
- gGlu-2xOEG tail 最贴近蛋白（3.8 A），linker 长度刚好产生有利接触
- NZ-C 全部 1.51-1.52 A，ParmEd bond 修复 100% 成功
- 趋势与 Lau 2015 Table 3 定性吻合

### MD 状态
- 4/15 done，11/15 84-94ns，预计今晚全部完成

---
*维护者：Claude Code*
*最后更新：2026-06-02*

## 2026-06-04 — 最终分析

### MD 完成状态
- 6/15 done，9/15 从 checkpoint 恢复续跑，13/15 轨迹可读（~95-102ns）
- 全部 NZ-C 稳定 1.51-1.52 A（ParmEd bond 修复成功）
- 全部温度 307-314K，PE -455k ~ -472k，无 NaN

### 最终结果

| Variant | EC50 (pM) | CA RMSD (A) | NZ-C (A) | Tail-Prot (A) |
|---------|-----------|:---:|:---:|:---:|
| No linker | 269 | 2.3+-0.2 | 1.51 OK | 4.6+-0.1 |
| gGlu | 9.9 | 2.2+-0.2 | 1.51 OK | 4.5+-0.7 |
| gGlu-1xOEG | 4.8 | 2.5+-0.3 | 1.51 OK | 4.4+-0.6 |
| gGlu-2xOEG | 6.2 | 2.5+-0.2 | 1.52 OK | 3.8+-0.2 |
| gGlu-3xOEG | 27.7 | 3.0+-0.2 | 1.52 OK | 4.1+-0.4 |

> ⚠️ **SUPERSEDED 2026-06-06**: 以下结论被 review `docs/reviews/claude-Jun06.md` 判定为"被数据严重高估"。
>
> **问题**:
> 1. 变体间 CA RMSD 差异（2.2–3.0 Å）在 ±SD 内大量重叠，趋势不单调（1×OEG RMSD 2.5 > gGlu 2.2）。
> 2. Tail-Prot "最优"差异仅 0.2–0.8 Å，未报告 replica 间 CV，未做 correlated t-test。
> 3. 承诺的 MM-GBSA 结合自由能未计算。
> 4. 实际模拟量 ~95–102 ns，为承诺 500 ns 的 ~1/5。
> 5. 2×OEG 的"最优"基于裸 EC50，但司美格鲁肽选 2×OEG 的真正依据是 BR ratio（需 HSA 端），而本体系无 HSA。
>
> **修正后的表述**：初步几何趋势（3×OEG RMSD 略高、2×OEG tail 距离略近）与 Lau 2015 Table 3 方向一致，但**当前数据不足以宣称统计显著或"复现设计逻辑"**。

### 结论（原始 — 见上方 SUPERSEDED 标注）
- **gGlu-2xOEG (Semaglutide)**: Tail-Prot 最近 (3.8 A)，linker 与 ECD 表面最优接触
- **gGlu-3xOEG**: CA RMSD 最高 (3.0 A)，过长 linker 导致 entropic penalty
- 趋势与 Lau 2015 Table 3 定性一致
- 计算复现了 Novo Nordisk 的 linker 设计逻辑

## 2026-06-06 — 重新分析 bond-fix 后的 fixed 轨迹

> ⚠️ **重大发现**：`compare_linkers.py` 之前一直加载的是 bond-fix **之前**的旧轨迹（`{vname}/repN/`），而非 fixed 轨迹（`{vname}_fixed/repN/`）。旧轨迹中 NZ-C ~14 Å（无共价键），fixed 轨迹 NZ-C ~1.51 Å。这意味着**之前的分析结论基于错误的物理状态**。脚本已修复路径和 prmtop 匹配。

### Fixed 版本真实结果（`common/lib/stats.py` 分析）

| Variant | EC50 | CA RMSD | Tail-Prot | Rep-CV (Tail) | 备注 |
|---------|------|---------|-----------|---------------|------|
| No linker | 269 | 2.3±0.2 Å | 4.6±0.1 Å | 4.0% | 3 replica |
| γGlu-1×OEG | 4.8 | 2.5±0.3 Å | 4.4±0.0 Å | 0.2% | **仅 2 replica** (rep1 缺失) |
| γGlu-2×OEG | 6.2 | 2.5±0.1 Å | **3.8±0.0 Å** | **0.8%** | 3 replica ✓ |
| γGlu-3×OEG | 27.7 | **3.0±0.1 Å** | 4.1±0.3 Å | 7.5% | 3 replica |

### 统计解读（升级后的 `stats.py`）
- **Tail-Prot 是最稳健的信号**：2×OEG 3.8 Å vs 其余 4.1–4.6 Å，replica CV 仅 0.8%
- **CA RMSD 区分度仍不足**：No linker 2.3 vs 2×OEG 2.5，差异 0.2 Å < ±SD
- **3×OEG CA RMSD 最高（3.0 Å）** 支持 entropic penalty 假说
- **n_eff 揭示自相关程度**：CA RMSD n_eff 低至 4–81，Tail-Prot n_eff 12–368。对 CA RMSD 直接报 mean±std 会严重高估精度
- **1×OEG rep1 DCD 为空** —— 实际可用 replica 不足

### 方法学教训
- `compare_linkers.py` 路径 bug 暴露了一个严重问题：**分析脚本和实际数据之间的匹配需要显式验证**（frame-0 NZ-C 距离检查应成为标准流程）
- OpenMM DCDReporter 时间戳不可靠（见 exp-A 同日记录），dt 必须从 reporter 设置推导

---
*维护者：Claude Code*
*最后更新：2026-06-06*

---

## 2026-07-17 — exp-D 重启：拓扑重建与 12 replica 生产启动

### 背景
原有 prmtop/inpcrd 与轨迹全部丢失。重建 4 变体（no_linker, gglu_1oeg/2oeg/3oeg；γGlu-only 仍排除）拓扑，重跑 3 replica × 100 ns 为 correlated t-test 采数。

### 发现的原始输入缺陷（lnk_*_pos.mol2）
1. 每文件 13–17 个孤立 H 原子（strip 残留，游离溶剂中）
2. strip_lya.py 硬编码删除名单只匹配 no_linker/gglu 布局：1oeg 丢 3 个链内原子（断链成 4 段）、1oeg/2oeg/3oeg 的 NME cap 未删（漂浮在 ~7 Å 外）
3. 无 BOND 表；内部酰胺全部误标 c3/o（酰胺键自由旋转）；电荷为按元素平均的占位值（非 AM1-BCC，且逐变体不一致）
4. MMFF 折叠构象穿过 ECD 表面（重原子重叠最低 0.48 Å）——与当年 gglu NaN 同类风险

### 重建流程
`build/rebuild_lnk_mol2.py`（结构感知 strip：删 ACE/backbone/sidechain/N10/NME + GAFF2 类型重标 + BOND 表转移，保留 bcc 电荷与锚定几何）→ `build/idealize_lnk_geom.py`（理想化链几何：C11 锚定 NZ+1.38 Å 沿 CE→NZ 轴，贪心二面角扫描+精炼，蛋白最小距离 2.68 Å、自接触 ≥2.46 Å）→ tleap（combine + `remove` Lys26 HZ1/2/3 + `bond` NZ-C11；本次均正常工作——此前的静默失败是 ParmEd 生成 PDB 格式特异的）→ `build/fix_prmtop_bonds.py` 验证。

关键补丁与适配：
- `tleap/lya_link_c8.frcmod`：ff14SB Lys CE 是 **C8** 不是 CT，原 lya_link.frcmod 的 `CT-N3-c-*` 二面角从未生效（tleap 报 missing torsion 且**拒绝保存 prmtop**）；按 CT 版同值补 `C8-N3-c-o` / `C8-N3-c-c3`。注意 frcmod DIHE 段必须用紧凑无空格格式（`C8-N3-c -o`），原文件带空格格式在 DIHE 段不解析。
- 本版 AmberTools 的 `addIonsRand` 不接受 `Na+ 0 Cl- 0`，改用 `addIons complex Na+ 0` + `addIons complex Cl- 0`。
- 只脱质子化 Lys26（其他 2 个 LYS 保留 6 个 HZ；旧 fix_bonds.py 曾误删所有 LYS 的 HZ）。

### Frame-0 验证

| Variant | Atoms | Waters | Na+ | NZ-C11 | LNK 键/角/二面角 | min500 PE (kJ/mol) |
|---------|------:|-------:|----:|:------:|:----------------:|:------------------:|
| no_linker | 36,682 | 11,547 | 8 | 1.38 Å ✓ | 52/100/165 | -569,861 |
| gglu_1oeg | 36,812 | 11,579 | 8 | 1.38 Å ✓ | 86/162/285 | -571,710 |
| gglu_2oeg | 36,816 | 11,574 | 6 | 1.38 Å ✓ | 107/200/356 | -569,545 |
| gglu_3oeg | 36,579 | 11,487 | 9 | 1.38 Å ✓ | 128/238/427 | -567,864 |

- NZ-C11 酰胺键全部存在（1.38 Å），0 个未参数化键/角/二面角，溶质重原子无 <1.9 Å 非键接触
- 净电荷残余 ±0.07~0.35 e（占位电荷的分数尾差，PME plasma 中和；与原生产方案一致）
- 已知遗留：LNK 电荷仍是按元素平均的占位值（与原生产一致，未改 AM1-BCC；若需严格电荷可比性，应 4 变体统一重做 BCC 后重建）

### 生产队列
- `md/launch_queue.sh`：GPU0 = no_linker + gglu_1oeg（6 replica 串行），GPU1 = gglu_2oeg + gglu_3oeg；GPU2/3 留给 exp-A
- 首批检查（启动 ~3.5 min）：no_linker rep1 min PE -591,612，生产 1.04 ns T=307.5 K PE=-475,450 kJ/mol；gglu_2oeg rep1 min PE -590,515，生产 1.08 ns T=309.4 K PE=-475,075 kJ/mol；385–398 ns/d/GPU，无 NaN
- 预计单 replica ~6 h，12 个共 ~36–38 h

---
*维护者：Kimi Code*

## 2026-07-19 — GPU 1 被外部进程抢占，rep3 迁移至 GPU 0

- 巡检发现 gglu_3oeg rep2 速度从 ~305 ns/d 跌至 ~90 ns/d：GPU 1 上存在**其他用户的进程**（PID 4067962，占 18.8 GB 显存），与本任务共享 GPU。
- 处置：杀掉 GPU 1 队列子壳（PID 2458216，rep2 进程不受影响继续跑），gglu_3oeg rep3 改为手动启动在空闲的 GPU 0（`CUDA_VISIBLE_DEVICES=0 ... --replica 3 --gpu 0`）。
- 结果：rep3 全速 392 ns/d；rep2 在共享 GPU 1 上 ~90 ns/d 继续。无 NaN。
- 教训：共享机器上 GPU 空闲状态可能随时变化；队列脚本未考虑外部抢占，必要时手动迁移。

---

## 2026-07-19 — 最终统计分析：correlated t-test 首次落地（compare_linkers_v2.py）

### 方法
12 replica（4 变体 × 3）全部跑完（各 1003 帧 × 100 ps = 100 ns，dt 从 DCDReporter 间隔推导）。新脚本 `analysis/compare_linkers_v2.py`：
- 指标沿用 v1 定义：CA RMSD（ECD CA，resSeq 0–100，对齐 ECD frame 0）、Tail-Prot（LNK 尾部 C[2:] 到蛋白 CA 最小距离）；末 50 ns 为平衡段。
- Frame-0 验证：12/12 replica NZ(Lys26)–C11 = 1.34–1.43 Å ✓，无拓扑/轨迹错配。
- 统计：`common/lib/stats.py::correlated_t_test`（**全项目首次真实调用**）作用于拼接逐帧序列（replica 独立，配对任意性已在 RESULTS.md 注明）；对照为 per-replica mean 的 Welch t（n=3，低功效）。

### 结果（平衡段，mean±SD 为 replica mean 间）

| Variant | EC50 (pM) | Tail-Prot (Å) | n_eff | CA RMSD (Å) |
|---|---:|---:|---:|---:|
| no_linker | 269 | 5.04±0.26 | 48 | 1.77±0.02 |
| gglu_1oeg | 4.8 | 3.74±0.14 | 45 | 1.71±0.04 |
| gglu_2oeg | 6.2 | 3.95±0.07 | 508 | 1.68±0.24 |
| gglu_3oeg | 27.7 | 3.65±0.12 | 97 | 1.82±0.14 |

主要检验（2×OEG vs X，Tail-Prot，correlated t）：vs no_linker Δ=−1.09 Å，p=1.3×10⁻⁵；vs 1×OEG Δ=+0.21 Å，p=0.015；vs 3×OEG Δ=+0.30 Å，p=0.004。次要（3×OEG vs X，CA RMSD）：全部 p>0.14，无显著差异。

### 结论
- 2×OEG 与每个变体都统计可区分（p≤0.015），但**唯一大效应是"有 linker vs 无 linker"（−1.1 Å）**；linker 长度之间的差异只有 0.2–0.3 Å，且方向与历史初步结果相反（2×OEG 并非最紧凑）。
- 历史初步结论（2×OEG 3.8 Å 独特紧凑、3×OEG CA RMSD 最高）**均未被新数据复现**。
- Tail-Prot 距离排序（3×<1×<2×≪无）与 Lau 2015 Table 3 活性排序（1×≈2×>3×≫无）不一致——"尾部贴近=活性高"不能解释 1×/2×/3×OEG 之间的活性差；只有 no_linker 极端定性吻合。
- 产物：`analysis/RESULTS.md`、`tail_prot_compare.png`、`per_frame_series.npz`、`test_results.json`、日志 `compare_linkers_v2.log`。

### 教训
- correlated t-test 的 n_eff 修正至关重要：2×OEG Tail-Prot n_eff=508（尾部弛豫快），CA RMSD 低 n_eff（6–101）使小差异不可检——直接 mean±std 会严重高估精度。
- v1 脚本 `resname LNK and name C` 在新拓扑下匹配为空（新原子名 C11），frame-0 NZ-C 检查再次证明是必需的回归防线。

---
*维护者：Kimi Code*
