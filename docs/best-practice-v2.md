# 司美格鲁肽 in silico 项目：最佳实践 v2

> 基于本项目（2026-05~06）的特有经验教训，补充 `docs/best-practice.md` 中未覆盖的内容。
> 来源：exp-A~F 实验日志、build 日志、review 文件、以及生产 MD 中的真实错误。
>
> **原则**：只记录本项目发生过的、有具体场景支撑的教训，不搬运通用最佳实践。

---

## 一、tleap / 拓扑构建

### 1. tleap 内部编号 ≠ 最终 prmtop 编号

- `bond`、`remove` 命令使用的是 tleap **加载顺序中的内部编号**，不是最终 PDB/prmtop 中的残基号。
- **场景**：`bond SYS.145.NZ SYS.157.C` 在日志中显示成功，但用户误以为是 atom 145 与 atom 157 成键，实际对应的是最终 prmtop 中的 resid 117 与 resid 129。
- **正确做法**：构建完成后，用 `parmed` 或 MDAnalysis 检查最终拓扑中的键列表，确认 `NZ` 与 `C` 之间确实存在 BondType。
- **代码**：
  ```python
  import parmed as pmd
  prm = pmd.load_file('system.prmtop')
  for bond in prm.bonds:
      if set([bond.atom1.name, bond.atom2.name]) == {'NZ', 'C'}:
          print(f"Bond: {bond.atom1.residue.name}{bond.atom1.residue.number}.{bond.atom1.name} -- "
                f"{bond.atom2.residue.name}{bond.atom2.residue.number}.{bond.atom2.name}")
  ```

### 2. tleap `bond` / `remove` 在某些 PDB 格式下静默失败

- 手动 Python f-string 生成的 PDB 可能导致 tleap 的 atom mask 解析失败，bond/remove 命令被静默忽略。
- **场景**：exp-D 中 `bond sys.117.NZ sys.129.C` 失败，导致 LNK 与肽之间无共价键（NZ-C 距离 14 Å），所有 linker 变体差异被掩盖。
- **正确做法**：
  1. 只用 **ParmEd `.save()` 生成的 PDB** 作为 tleap 输入（其对齐格式被 tleap mask 引擎正确解析）。
  2. bond 命令执行后，**立即检查最终 prmtop 的键列表**。
  3. 如果 tleap bond 不可靠，改用 **ParmEd `BondType(k, req)` 在 prmtop 层面直接添加键**。

### 3. 跨力场残基连接需要显式 cross-force-field 参数

- 当 ff14SB 残基（如 LYS）与 GAFF2 小分子（如 LNK）形成共价键时，tleap 不会自动识别 N3(amber)-c(gaff2) 这类跨力场键参数。
- **场景**：LYS NZ — LNK C 的酰胺键需要显式 `BondType(427.0, 1.38)`，否则模拟中键长会漂移到 14 Å。
- **正确做法**：在 frcmod 中补充跨力场的 bond/angle/dihedral 参数，或在 ParmEd 中直接注入 BondType。

### 4. tleap 跳过 H 原子可避免 NGLY/CGLY/ACE 模板冲突

- 当 N-terminal 残基同时携带 H 原子和 ACE cap 时，tleap 的 NGLY 模板坐标不匹配导致 FATAL 错误。
- **场景**：exp-F ECD v2 构建中，HSA + 肽 + ECD 组合时 N-terminal GLY 报错。
- **正确做法**：在 build 脚本中 `if element == 'H': continue`，让 tleap 自动补全 H 原子。

### 5. `loadAmberParams` 后必须加载对应力场源文件

- 加载 GAFF2 原子类型的 mol2 前，必须 `source leaprc.gaff2`，否则 `c`, `c3`, `os` 等类型找不到 vdW 参数。
- **场景**：exp-F Phase 0 中遗漏 `leaprc.gaff2`，LNK 原子因无 vdW 参数导致 NaN。

---

## 二、mol2 / 小分子参数化

### 6. 长链脂肪酸 mol2 容易混入孤立原子

- C18 二酸链长，mol2 中原子编号容易出错，末端几个原子最容易遗漏或重复。
- **场景**：`lnk_noh_fixed.mol2` 含 4 个与主链断开的原子（C58/O59/N60/C61），导致 72% 帧 PE 跳变 >1000 kJ/mol。
- **正确做法**（三重检查）：
  1. **原子数核对**：与 SMILES/预期结构对比（LNK 应为 43 原子，不是 47）。
  2. **末端连通性**：用 parmed 检查末端碳（C55）的键列表，确认有 C54-C55 和 C55-C56（或 C55-O）。
  3. **可视化**：VMD 加载后旋转查看链末端，孤立原子在 ~40 Å 外一眼可见。
- **教训**：能量异常（PE 跳变 >1000 kJ/mol）是拓扑错误的第一信号，应立即停跑检查。

### 7. antechamber/sqm 在非标准分子上不可靠

- sqm 对 γGlu-2×OEG-C18 这类复杂分子可能产生断裂几何：CA-CB 键长 2.61 Å（标准 1.54 Å），热函 +257 kcal/mol。
- **场景**：exp-C 中手工 GAFF2 电荷（-0.80~+0.70）过于极化，导致 NaN；改用 AM1-BCC 后（-0.35~+0.10）稳定运行。
- **正确做法**：
  1. 对复杂非标准分子，优先用 AM1-BCC（`antechamber -nc 0 -c bcc`）。
  2. sqm 输出必须检查：最终热函应 < +50 kcal/mol，键长在标准范围内。
  3. 如果 antechamber 失败（sqm/bondtype 子程序不可用），考虑 RDKit + Open Babel 替代管线。

### 8. 手工电荷必须检查总电荷和极化范围

- 手工分配的 GAFF2 电荷范围 -0.80~+0.70 明显过极化，AM1-BCC 的 -0.35~+0.10 更合理。
- **场景**：全局 charge scaling（+6.75 → -1）导致符号反转，NaN at 0.3 ns。
- **正确做法**：
  - 总电荷必须精确到 0.0001（`sum(charges)`）。
  - 单个原子电荷绝对值不应超过 ~0.5（除非是离子）。
  - 偶数电子体系确保 SQM 收敛。

---

## 三、膜体系构建

### 9. packmol-memgen `renumber=True` 导致残基号冲突

- `MembraneParams.pdb_reindex(renumber=True)` 将蛋白残基按链重新编号从 1 开始。
- `charmmlipid2amber.py` 仅以 `(chain, resnum)` 识别残基，导致脂质与蛋白残基合并、肽键断裂。
- **场景**：exp-F 膜体系构建中，packmol-memgen 生成的 PDB 无法直接用于 tleap。
- **正确做法**：使用自定义 merge pipeline（Python + scipy.cKDTree 移除重叠脂质），绕过 packmol-memgen 的 reindex。

### 10. PDB 固定宽度格式：4-char resname 导致列偏移

- Python f-string `{resname:>3s}{chain:>1s}` 在 resname 为 4 字符（如 `WATA`）时，chainID 移至 col 21，坐标列整体左偏 1 列。
- **场景**：相邻负坐标合并为 `-43.108-107.555`。
- **正确做法**：resname 和 chainID 之间保留**显式空格**：`{resname:>3s} {chain:>1s}`。

### 11. 膜体系规模以实际构建结果为准（本项目膜体系 312k atoms）

- 直觉认为"加膜后原子数更多"，但实际情况取决于 box 尺寸。
- **场景**：溶剂化体系 363k atoms（oct box 171.6 Å），膜体系 **312,476 atoms**（orthorhombic box 更小）。
- **教训**：不要预设膜体系一定更大，以实际构建结果为准。后续性能估算必须基于真实原子数。

### 12. 膜环境是跨膜受体稳定的必要条件

- 无膜溶剂化体系中，GLP-1R TMD 在 25 ns 内漂移 32 Å（RMSD），完全丧失结构意义。
- 膜体系中受体 RMSD 仅 3 Å，肽段 <1.2 Å。
- **教训**：分析跨膜蛋白时，"无膜"数据直接作废，无需犹豫。

---

## 四、轨迹与数据完整性

### 13. OpenMM DCD 时间戳不可靠

- `DCDReporter` 写入的时间戳单位与 mdtraj 解析期望不一致（fs→ps 转换问题），dt 被低估 ~100 倍。
- **场景**：mdtraj 读取显示 dt=1.0 ps，实际 reporter 设置是 50000 steps × 2 fs = 100 ps。
- **正确做法**：分析脚本中 `dt_ns` 从 reporter 间隔**硬编码**推导，绝不从 DCD 时间戳读取。

### 14. OpenMM 正在写入的 DCD 可被 MDAnalysis 安全读取

- 生产 MD 期间可以实时分析 DCD，但帧数会随写入进度动态增加。
- **正确做法**：多个 DCD 并行分析时，取各 universe 的 `len(trajectory)` 最小值，避免维度不匹配。

### 15. DCD 文件在 MD 运行中读取可能损坏

- 如果多个进程同时访问同一 DCD（如重复启动 launcher），文件头或帧体可能损坏。
- **场景**：exp-C rep1 DCD 因前期重复进程损坏，分析脚本需容错跳过。
- **正确做法**：分析脚本使用 `md.iterload`（mdtraj）或 `Universe()` + `try/except`（MDAnalysis），损坏帧 graceful 跳过。

### 16. XML 缓存有 bug（OpenMM System 反序列化后 addForce hang）

- 反序列化的 System 对象在添加 CustomExternalForce 时，内部 reindex 可能死锁（CPU 0%, GPU 100% 无进展）。
- **场景**：exp-C 中 XML 反序列化（31 MB）后 addForce  hang。
- **正确做法**：从 prmtop 重建系统（`createSystem` 仅 3s）比 XML 缓存（40+s 加载 + 风险）更可靠。

### 17. GROMACS `md.xtc` 默认保存 wrapped 坐标，直接读入 MDAnalysis 会导致 RMSD 严重低估

- GROMACS 压缩轨迹默认保存 PBC-wrapped 坐标，跨 box 边界的原子会被折回，使柔性区域（如 linker、loop）看起来"没有移动"。
- **场景**：exp-F GROMACS 膜体系中，LNK RMSD 直接从 `md.xtc` 读取仅 3.9 Å，修正 PBC 后实际为 **10.9 Å**——低估了 **~2.8×**。
- **正确做法**：
  1. 分析前先用 `gmx trjconv -pbc whole` 生成完整分子轨迹：
     ```bash
     echo "0" | gmx trjconv -s md.tpr -f md.xtc -o md_whole.xtc -pbc whole -n system.ndx
     ```
  2. 以 `md_whole.xtc` 作为 MDAnalysis/mdtraj 的输入。
  3. 对膜体系，`whole` 比 `nojump` 更适合，因为 `nojump` 会让整个系统随 COM 漂移出 box。
- **教训**：分析 GROMACS 轨迹时，**永远不要直接信任原始 xtc 的坐标**，PBC 处理必须是标准流程的第一步。

---

## 五、MD 引擎性能与 GPU 分配

### 27. 大型膜系统（312k atoms）GROMACS 2026 比 OpenMM 快约 2.3×

- **实测数据**（RTX 3090，2 fs 步长，各自独占 GPU）：
  - GROMACS 2026：~~81–100~~ **96.4 ns/day**（`-nb gpu -pme gpu -update gpu`，PP/PME/约束全 GPU）
  - OpenMM 7.7+：~~约 95~~ **42.1 ns/day**（CUDA mixed precision，Monte Carlo barostat）
- **关键发现**：系统实际为 **312,476 atoms**，此前误记为 ~140k。基于 312k 的 2.3× 差距是合理的（GROMACS 的 Verlet 列表和 PME 优化在大系统上优势放大）。
- **教训**：
  1. 性能估算前必须确认**真实原子数**。
  2. GROMACS 对 300k+ 原子膜系统的优势可达 2× 以上，跨引擎比较时不能假设同速。

### 28. 多 GPU 环境下必须显式分配设备，避免隐性竞争

- **场景**：4× RTX 3090 机器上，OpenMM 默认取 device 0，GROMACS 默认也取 device 0，若同时启动会互相抢占导致双方速度暴跌。
- **正确做法**：
  - OpenMM：`platform_props = {"CudaDeviceIndex": "1"}`
  - GROMACS：`export CUDA_VISIBLE_DEVICES=2` + `gmx mdrun -nb gpu -pme gpu -update gpu`
- **验证**：`nvidia-smi` 确认两个进程分别落在不同 GPU，utilization 各自 ~99%，无竞争。

### 29. 速度必须用 wall-time 独立验证，不能只看引擎自报

- OpenMM 的 `StateDataReporter speed` 与独立 wall-time 计算（`ns_run / hours_wall × 24`）吻合（42.3 vs 42.1 ns/day）。
- GROMACS 没有实时 speed 列，必须用 `tpr` 创建时间到 `log` 修改时间推算；早期曾误将 `ns/h` 当成 `ns/day` 导致低估。
- **教训**：
  - 对 GROMACS，以 `md.log`/`npt.log` 的 `Modify` 时间减去 `md.tpr`/`npt.tpr` 的 `Birth` 时间作为 wall-time。
  - 计算时先统一为 `ns/h`，再**显式 ×24** 得到 `ns/day`。

### 30. OpenMM 的 Python checkpoint 每 1 ns 写入有一定 overhead

- **场景**：OpenMM 每 1 ns 调用 `context.createCheckpoint()` 写入 46 MB 文件，在 312k 原子系统上单次写入耗时数秒。
- **教训**：长生产 run 的 checkpoint 间隔可以适当放宽（如 5 ns），减少 I/O 阻塞；GROMACS 的 `.cpt` 由引擎内部高效管理，无需担心。

---

## 六、分析方法论

### 17. 柔性 linker 体系必须区分内部稳定性与域间运动

- `superposition=False` 的 RMSD 会将正常的域间平移/旋转（几十 Å）误判为"结构破坏"。
- **场景**：ECD v2 中 superposition=False RMSD 高达 35-55 Å，改用 superposition=True 后内部 RMSD 仅 1-3 Å。
- **正确做法**：
  - **内部稳定性**：`superposition=True` RMSD（结构是否保持折叠）。
  - **域间运动**：COM 距离（各组分质心之间的距离变化）。
  - **Linker 柔性**：Rg（半径回转）或端到端距离。
  - 三个指标必须同时报告，缺一不可。

### 18. RMSF 必须采用 component-wise alignment

- 多组分复合物（如 ECD+Peptide+HSA）中，全局对齐会将域间运动 artifactually 计入 RMSF。
- **场景**：ECD v2 中联合对齐导致 RMSF 被高估数倍。
- **正确做法**：每个组分独立对齐后分别计算 RMSF：
  ```python
  for sel in ['resid 1-100', 'resid 101-128', 'resid 130-711']:
      u2 = mda.Universe(PRMTOP, DCD)
      align.AlignTraj(u2, ref, select=f'name CA and {sel}', in_memory=True).run(step=STEP)
      rmsf = RMSF(u2.select_atoms(f'name CA and {sel}')).run(step=STEP)
  ```

### 19. MM-GBSA 需要干态轨迹

- MMPBSA.py（cpptraj 后端）要求轨迹原子数与干态 prmtop 匹配。
- **场景**：直接传入溶剂化 DCD（144814 原子）到干态 prmtop（2065 原子）导致 cpptraj 报错。
- **正确做法**：先用 cpptraj 提取干态轨迹：
  ```
  parm system.prmtop
  trajin prod.dcd
  strip :WAT,Na+,130-711  # 去溶剂 + 去掉不需要的组分
  trajout dry.dcd dcd
  run
  ```

### 20. ante-MMPBSA.py 的 Amber mask 从 1 开始计数

- parmed 的 residue 索引从 0 开始，但 Amber mask `:N` 从 1 开始。
- **场景**：`:1-100` 只选了前 100 个残基（THR0-GLU99），漏掉了 GLY100，导致配体包含 HSA 的大部分。
- **正确做法**：
  1. 先用 parmed 打印残基列表，确认目标残基的 **Amber mask 编号**（= parmed 索引 + 1）。
  2. 拆分后用 parmed 验证：检查受体和配体的原子数、首/末残基名是否与预期一致。

---

## 七、晶体结构与模型构建

### 21. 所有 GLP-1R ECD-肽晶体结构中 ECD 与肽均分离 25-40 Å

- 这不是晶格堆积造成，而是活性态 GLP-1R 的真实结构特征（ECD 远离 TMD 正位点）。
- **场景**：exp-B 中 3IOL/4ZGM/6GB1/7KI0 全部验证，min dist = 24.9-38.6 Å。
- **教训**：不要试图用对接/refinement 桥接这个间隙——FlexPepDock 对 >10 残基肽不可靠。

### 22. 游离脂肪酸不能替代 linker-FA 复合体

- linker 接上后，FA 从 HSA FA3 口袋逃逸（32-41 Å），无论末端电荷如何。
- **场景**：exp-C 中游离 C18 二酸锚定 ARG482（2.8 Å），但 linker-C18 逃逸至 32 Å；ACE 封端（消除 NH₃⁺）后仍逃逸至 38-41 Å。
- **物理解释**：OEG 单元的亲水性无法容忍 FA3 疏水口袋的埋藏环境。
- **教训**：参数化时必须测试完整 linker-FA 复合体，不能仅用游离 FA 推断行为。

---

## 八、文档与项目管理

### 23. 文档声明不能跑在实际状态前面

- **场景 1**：exp-D commit message 写 "qualitative match to Lau 2015"，实际差异在误差棒内、无统计检验、无 MM-GBSA。
- **场景 2**：Jun01-response 声称 barostat "已修"，实际 `common/scripts/run_md.py` 未改。
- **场景 3**：roadmap 标 exp-D "完成"，但 "待启动"一行仍写 "exp-D"。
- **正确做法**：
  - 任何结论声明必须附带置信度标签（"preliminary / within error / not statistically significant"）。
  - "已修"必须有代码/文件的 diff 证据，不能只有口头承诺。
  - 被推翻的旧结论在原条目顶部标注 `> ⚠️ SUPERSEDED YYYY-MM-DD: 见XX修正`。

### 24. tasks.md 不更新会变成僵尸清单

- **场景**：exp-A tasks.md 全部是 ⬜ 未做，而 exp-log 显示 200 ns 生产 + 全量分析早已完成。
- **正确做法**：要么定期同步 tasks.md，要么废弃 tasks.md 改用 exp-log.md 作为唯一进度源。

### 25. 分析脚本路径错误可能导致结论基于错误物理状态

- **场景**：`compare_linkers.py` 加载了 bond-fix 之前的旧轨迹（NZ-C ~14 Å），而非 fixed 版本（NZ-C ~1.51 Å），得出"变体间无差异"的错误结论。
- **正确做法**：frame-0 验证（关键距离检查）必须成为标准流程，确认加载的是正确的拓扑和轨迹。

### 26. 生产 MD 完成后必须立即归档原始拓扑，不能在同目录重建

- **场景**：exp-F 膜体系第一次生产 MD（110 ns）完成后，用户在 `membrane_build/` 目录重新执行 tleap build 流程。新的 `system_final.prmtop` 直接覆盖了旧版本，导致旧 prmtop（312,501 atoms）永久丢失。
- **后果**：
  1. 旧 checkpoint（55M steps）虽然能加载，但坐标已损坏（PE = 67.6 亿 kJ/mol，原子散布 2800 Å），第一步即 NaN。
  2. 旧 DCD（68 ns）无法与任何现有 prmtop 匹配（原子顺序不同），PE 计算异常（65 亿 kJ/mol）。
  3. **110 ns 轨迹数据完全作废**，无法延续，也无法用正确能量重新分析。
- **根本原因**：旧 run 使用的 `membrane_build/system_final.prmtop` 是**唯一**能正确解释旧 DCD/checkpoint 的拓扑，没有备份即被覆盖。
- **正确做法**：
  1. **生产 run 启动前**，将构建好的 `system_final.{prmtop,inpcrd,pdb}` 复制到**归档目录**（如 `archive/system_v1/`）。
  2. **生产 run 完成后**，将 `prmtop + inpcrd + 所有 checkpoint + DCD + log` 打包到 `archive/run_v1/`。
  3. 只有当归档完成且验证通过后，才允许在同一目录重建系统。
  4. 版本命名规则：`system_v1.prmtop`、`system_v2.prmtop`，绝不在同一文件名上覆盖。

---

## 九、安全

### 26. API key 安全问题会以新文件名复发

- `cc-ds.sh` 被 gitignore 后，新出现的 `cc.sh` 仍含明文 token 且未被 gitignore。
- **根治方案**：所有启动脚本统一从 `.env`/环境变量读取，杜绝"换个文件名又泄露"。

---

## 十、快速检查清单（本项目特有）

### 构建后
- [ ] mol2 原子数 vs 预期 SMILES 一致
- [ ] mol2 末端连通性（parmed 键列表）
- [ ] tleap bond/remove 后在最终 prmtop 中验证键存在
- [ ] 跨力场连接处有显式 bond/angle 参数
- [ ] 电荷总和精确到 0.0001，单个原子 |charge| < 0.5
- [ ] VMD 可视化确认无孤立原子

### MD 启动前
- [ ] frame-0 关键距离检查（如 NZ-C 酰胺键、催化距离、盐桥）
- [ ] 初始 PE 合理（无 NaN/10²¹）
- [ ] dt 硬编码在分析脚本中，不从 DCD 时间戳推导

### 分析前
- [ ] 确认加载的是正确的拓扑和轨迹（非旧版本）
- [ ] PBC 已处理（GROMACS 轨迹）
- [ ] 柔性 linker 体系同时报告：内部 RMSD + COM 距离 + Rg
- [ ] RMSF 采用 component-wise alignment
- [ ] MM-GBSA 使用干态轨迹

---

*最后更新：2026-06-09（GROMACS 96.4 vs OpenMM 42.1 ns/day 实测数据已补充）*
*来源：semaglutide-in-silico 项目 exp-A~F 实验日志、构建日志、review 文件*
