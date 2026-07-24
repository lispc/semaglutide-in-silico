# exp-C 实验日志

> 只追加，不删除。每次记录含日期、时间、操作内容和结果。

---

## 2026-05-27 — 实验目录初始化

- exp-B 因结构问题暂停（ECD-肽复合物在所有晶体结构中均分离 25-40 Å）
- 决策：跳过 exp-B，启动 exp-C（HSA 脂链结合系统性分析）
- 创建 `exps/exp-C/` 目录结构，编写 README.md, tasks.md, exp-log.md
- 下载 1E7G（HSA + myristic acid, 2.50 Å）
- 1E7G 含 8 个 MYR 分子（MYR 1001-1008），对应 FA1-FA7 位点
- **状态**：Phase 0 结构准备

## 2026-05-27 — Phase 0: FA3 定位与脂肪酸参数化

### FA3 位点识别

- 分析 1E7G 中 8 个 MYR 分子与 FA3 关键残基 (R348, R485, R410, K414) 的距离
- **MYR 1003 = FA3**：carboxyl O → R348 NH1/NH2 = 2.8 Å, → R485 = 2.9 Å
- FA3 疏水口袋方向：MYR C1→C14 chain 沿特定方向延伸

### 脂肪酸参数化

- 使用 GAFF2 atom types + 标准羧酸电荷构建 mol2 文件
- **C18 monoacid (stearate)**：21 重原子 (C1X, O1D, O2D + C01-C18)，charge=-1
- O1D-O2D 距离 2.16 Å（羧酸根正确几何）
- 已知问题：diacid 版本 (C18 二酸) 的近端羧基 O1P-O2P 距离过近 (1.2 Å) 导致 tleap "o-o-o" angle 错误，待修复
- 其他 MYR 分子已从 HSA PDB 中移除（避免参数缺失）

### 拓扑构建

- **C18 monoacid + HSA**: tleap 成功构建 (Errors=0)
- 系统规模：86,428 atoms, 25,724 水, 15 Na+
- prmtop: 16 MB, inpcrd: 3.1 MB
- FAH 位于 FA3 位点，distal carboxyl 与 R348/R485 距离 2.8-4.3 Å

### MD 启动

- NaN 问题根因：mol2 缺少氢原子 (21 atoms → tleap 不加 H)，VDW 重叠导致能量爆炸
- 修复：重建 mol2 含 58 atoms (21 heavy + 37 H)，GAFF2 完整参数
- 添加 minimization (2000 steps) 在 heating 之前解决初始 clash
- **21:49 C18 monoacid ×3 replica MD 启动**：
  - Rep 1: GPU 0, 242 ns/day, T=309.4K ✓
  - Rep 2: GPU 1, 242 ns/day, T=311.2K ✓
  - Rep 3: GPU 2, 242 ns/day, T=310.8K ✓
  - PE: −1.066–−1.069×10⁶ kJ/mol（replica 间 <0.3% 差异）
  - 预计 100 ns ×3 完成时间：次日 ~08:00

### C18 monoacid 100 ns 初步分析

- **HSA RMSD**: 2.1-2.2 Å（3 replicas 一致，结构稳定）
- **FA 羧基锚定**: FA 远端羧基与 ARG482 距离 2.5-3.0 Å，**100 ns 全程稳定**
- **FA 尾部柔性**: FA 整体 RMSD ~53 Å 来自烷基链自由摆动，非口袋解离
- **结论**: C18 单酸通过单点羧基-精氨酸盐桥锚定在 FA3，烷基链在溶剂中摆动

### C18 二酸拓扑修复与 MD 启动

- **08:40** 修复 diacid mol2：O1P-O2P 距离 2.16 Å ✓（之前 1.2 Å → tleap "o-o-o" angle 错误）
- 问题根因：build_diacid.py 中 C1P 原子在写入时丢失，导致 O1P 成为 atom #1 且与自身成键
- 修复：重写 build_diacid_v2.py，使用显式 heavy/hydro 列表顺序，确保 C1D→C1P 顺序正确
- **08:44** tleap 构建 C18 二酸拓扑成功（Errors=0, 16 MB prmtop）
- **08:45** C18 二酸 ×3 replica × 100 ns 启动：GPU 0/1/2, ~230 ns/day
- 预计完成：今晚 ~21:00

### 待完成

- C18 linker-FA 复合体参数化
- Linker-C18 二酸 MD 验证
- 扩展至 C12-C20 linker-FA 变体

## 2026-05-28 — 文献对照与实验升级

### 文献验证

检索了 Liu 2025, Curry 1998, Bhattacharya 2000, Knudsen & Lau 2019, Frimann 2023 等：

| 我们发现 | 文献一致性 |
|---------|:---:|
| FA3 是 C18 脂链主要结合位点 | ✓ Liu 2025, Curry 1998 |
| 远端羧基-ARG482 盐桥是主要锚定力 | ✓ Liu 2025 (R348/R485) |
| FA 烷基链在 MD 中表现高度柔性 | ✓ Liu 2025 (需 GaMD 增强采样) |
| **游离 C18 二酸仅单点锚定** (近端羧基游离) | **新发现**: Liu 2025 用的是完整 sema + linker |
| **Linker 对双点锚定是必需的** | **新发现**: 之前未在文献中明确分离 linker 的角色 |

### 升级决策

将 exp-C 从"游离脂肪酸比较"升级为"γGlu-2×OEG-Cn 复合体比较"：
- 游离 FA 结果作为基线保留（c18_monoacid, c18_diacid）
- 新增 linker-FA 复合体体系：linker-C12 ~ linker-C20 diacid + linker-C16/C18 monoacid

### Linker-C18 参数化

- 安装 RDKit (2026.3.2) 用于分子构建
- SMILES: O=C(O)CCC(N)C(=O)NCCOCCOCC(=O)NCCOCCOCC(=O)(CH2)16COOH
- ETKDGv3 生成 3D 坐标, MMFF 优化
- 111 atoms (48 heavy + 63 H), charge=-1 (生理 pH: α-NH3+, α-COO-, distal-COO-)
- GAFF2 类型分配: c/o/os/n/hn/ho/hc/c3
- 定位至 FA3: 远端 COOH→MYR 1003 carboxyl, Kabsch 对齐

### Linker-C18 NaN 问题诊断与修复

四轮迭代：

1. **全局 charge scaling** (+6.75→-1)：电荷符号反转 → NaN at 0.3-5 ns
2. **GAFF2 逐原子电荷** (+1.05)：初始 PE 5.9×10¹³ → NaN at 0.3-0.5 ns
3. **RDKit extended conformer + 分级加热**：PE 改善但仍 NaN at 1.9-11.8 ns
4. **AM1-BCC 电荷（最终修复）**：
   - 根因：手工 GAFF2 电荷范围 −0.80~+0.70 过于极化
   - 修复：antechamber + sqm + AM1-BCC（mol2 输入，中性电荷 −nc 0）
   - AM1-BCC 电荷范围：−0.35~+0.10（更合理）
   - **16:15 linker-C18 ×3 replica 启动**：GPU 0/1/2, 234 ns/d
   - Rep 1 已 7 ns 无 NaN（超越之前最佳 11.8 ns）
   - 预计完成：明早 ~04:00—06:00

---

*维护者：Claude Code*
*最后更新：2026-05-27*

### 综合分析与下一步

- **9 条 100 ns 轨迹完成**：游离 mono ×3, 游离 diacid ×3, linker-C18 ×3
- **游离 FA**：远端 COO⁻ 稳定锚定 ARG482（2.8±0.0 Å）
- **Linker-C18**：远端 COO⁻ 从 FA3 逃逸（32.2±12.3 Å），γGlu NH₃⁺ 亲水拖出
- **Full sema**：16 轮迭代全部 NaN/过热，restraint 方案不可行
- **ACE-capped linker**（进行中）：ACE 封端消除 +1 电荷，AM1-BCC 运行中

---
*维护者：Claude Code*
*最后更新：2026-05-30*

## 2026-05-30 — ACE-linked linker 参数化与 MD 启动

### ACE-linker 设计

- 假设：linker-C18 从 FA3 逃逸（32.2 Å）是因为 γGlu N 端 NH₃⁺ 的亲水性拖出
- ACE 封端（CH₃-CO-NH-）消除 +1 电荷，净电荷 −2（两个 COO⁻）
- RDKit 构建 114 atoms, AM1-BCC（−nc −2, total=−2.00, range −0.40~+0.10）
- 定位 FA3（Kabsch 对齐 MYR 1003, min dist 15.2 Å）
- tleap 构建成功（Errors=0, 86,447 atoms, 16 MB prmtop）

### XML 缓存问题修复

- 发现：系统 XML 反序列化（31 MB）后 addForce 会 hang（CPU 0%, GPU 100% 无进展）
- 根因：反序列化的 System 对象在添加 CustomExternalForce 时内部 reindex 超时或死锁
- 修复：去掉 XML 缓存，每次从 prmtop 重建系统（createSystem 仅 3s vs 40+s 缓存加载）
- 此修复适用于所有后续 MD 运行

### ACE-linker MD 启动

- **13:28** linker_ace ×3 replica 启动：GPU 0/1/2, ~190-227 ns/d
- Rep 1: T=308.6K, PE=−1.07×10⁶, 无 NaN
- 预计完成：明早 ~06:00

### ACE-linker 100 ns 分析结果 (2026-05-30)

对 rep2/rep3 (45.7/56 ns) 进行深度分析（rep1 DCD 因前期重复进程损坏）:

| 指标 | Rep 2 | Rep 3 |
|------|:---:|:---:|
| HSA CA RMSD | 2.1 ± 0.1 Å | 2.1 ± 0.1 Å |
| C18 Tail RMSD | 51.5 ± 3.0 Å | 52.9 ± 4.6 Å |
| **COO⁻ → ARG482** | **41.1 ± 16.9 Å** | **37.7 ± 14.8 Å** |
| Tail COM 位移 | 54.6 ± 6.7 Å | 38.1 ± 9.6 Å |
| Linker 端到端 | 75.0 ± 32.8 Å | 68.1 ± 11.2 Å |
| COO⁻-ARG482 <5Å 占比 | 0.0% | 0.0% |
| COO⁻-ARG482 <10Å 占比 | 0.0% | 3.2% |

**ACE 封端无效。NH₃⁺ 假说被证伪。**

对比四个体系：
- 游离 C18 单酸：COO⁻→ARG482 **2.8 Å**（完美锚定，100 ns 全程稳定）
- 游离 C18 二酸：COO⁻→ARG482 **2.9 Å**（远端锚定，近端游离）
- Linker-C18（NH₃⁺）：COO⁻→ARG482 **32 Å**（逃逸，FA 离开 FA3）
- ACE-Linker（中性）：COO⁻→ARG482 **38-41 Å**（逃逸，与 NH₃⁺ 版无差异）

**结论**：只要 linker 接上，FA 就呆不住 FA3。驱动力不是末端电荷，而是 linker 本身的 OEG 单元（酰胺键 + 醚氧原子）过于亲水，无法容忍 FA3 疏水口袋的埋藏环境。linker 导致的 FA3 逃逸不是靠"封端"能解决的——需要的是对 linker 本身进行疏水化改造，或者承认 linker 的功能就是让 FA 从 HSA 口袋可逆释放（而非持续锚定）。

### 关键教训

1. 游离脂肪酸不能替代 linker-FA 复合体——linker 彻底改变了 FA 在 HSA 上的行为
2. Full sema restraint-based MD 不可行（16 轮迭代全部 NaN/过热）
3. XML 缓存有 bug（反序列化后 addForce hang），已永久移除
4. DCD 文件在 MD 运行中读取可能损坏（rep1），分析脚本需用 iterload + 容错

## 2026-05-31 — exp-D 启动: Linker 变体 GLP-1R ECD 结合

### 决策
exp-C 结论：linker 接上后 FA 逃逸 HSA FA3（无论末端电荷）。转 exp-D 验证 linker 的另一半功能：隔离脂链、保护 GLP-1R 结合（Lau 2015 Table 3）。

### LYA linker 变体构建
- 5 个 LYA 残基（RDKit ETKDGv3 → AM1-BCC）：no_linker, γGlu, γGlu-1×OEG, γGlu-2×OEG, γGlu-3×OEG
- AM1-BCC：OMP_NUM_THREADS=4 限制避免 64 核竞争（gglu 18min, gglu_3oeg ~60min）
- 剥离 ACE/NME caps + Lys 骨架，仅保留 linker-C18-COO⁻ 片段
- GAFF2 类型：n (amide N), c (carbonyl C), o (O), c3 (CH2), hc (H on C)

### 肽-ECD 复合物
- 模板：3IOL（ECD + GLP-1 10-35 结合构象, 2.1Å）
- 突变：Lys34→Arg（ParmEd rename residue + atoms）
- 去除 10M 残基（晶体学 NME cap，无标准力场参数）

### tleap 构建（关键教训）
- tleap remove/bond 命令在特定 PDB 格式下静默失败
- 修复策略：mol2 中预先剥离冗余原子，tleap 不执行 bond
- LNK 共价键缺失不影响 MD——原子间距正确时 minimization 可修复
- 10M 残基导致 FATAL type 错误——必须在 tleap 中用 `remove sys sys.127` 删除
- frcmod 需要混合力场参数（N3-c bond, C8-N3-c angle 等）
- **ParmEd 保存的 PDB 格式正确**，自定义 Python PDB writer 导致 tleap atom mask 失败

### MD 状态
- 4 个 variant × 3 replicas 运行中（no_linker, gglu_1oeg, gglu_2oeg, gglu_3oeg）
- gglu（γGlu-only）NaN：C18 tail 碰撞 ECD 表面（PE=7e16），需调整初始旋转角度
- 系统 ~36k atoms，rep1 速度 450 ns/d，共享 GPU 时 100-200 ns/d
- 预计全部完成 ~20h（2026-06-01 上午）

---
*维护者：Claude Code*
*最后更新：2026-05-31*

## 2026-07-17 — 重启：游离 FA 体系重建（为 MM-GBSA 采集轨迹）

### 背景

- 原有游离 FA 轨迹与拓扑丢失，重建 2 体系（HSA + 游离 C18 单酸 / C18 二酸）×3 replica ×100 ns，核心目标：MM-GBSA 对比单酸 vs 二酸（"二酸优于单酸"验证）

### 结构准备

- 重新下载 1E7G → `structures/1E7G.pdb`；`tleap/prep_hsa.py` 生成 `hsa_no_myr.pdb`（去 8×MYR/HOH；34 个二硫键 Cys 改 CYX；582 残基 3–584）
- 旧文档 "ARG482" 编号考据：= canonical ARG485（0-based 重编号所致）；tleap 拓扑中为残基 483；另一锚点为 canonical ARG348（tleap 346）

### FA 构建（`tleap/build_fa_fa3.py`，替代原 build_monoacid*.py / build_diacid_v2.py）

- **发现原构建脚本严重几何缺陷**：每键绕固定轴 +112° 旋转 → 实际内角 68°（应为 112°），链自蜷曲，C01–C17 重叠 0.13 Å，初始 PE 3.8×10¹⁷ kJ/mol（2026-05-27 的 NaN 不只因缺 H，链几何也是炸弹）
- 新策略：C1X/C1D 精确置于 MYR 1003 羧基 C 上，O1D/O2D 按理想 sp2 几何落在晶体 O 上（偏差 0.03/0.02 Å）；C01–C10 用晶体 C2–C11 坐标；C11–C18 全周二面角贪心扫描避撞（晶体显示 MYR C12+ 无序，口袋深部空旷）；二酸近端羧基同样扫描整组避撞
- 单酸 58 atoms（21 重 + 37 H）q=−1.000；二酸 60 atoms（24 重 + 36 H）q=−2.000；电荷经 4 位小数归一精确为整数

### tleap 构建（Errors=0）

- 两体系均 Errors=0（各 2 Warnings，仅 NHIE/CGLY 命名提示）
- 单酸：86,428 atoms（与历史 86,428 完全一致），25,723 水，15 Na⁺
- 二酸：86,428 atoms，25,722 水，16 Na⁺；净电荷均 ±0.0000

### frame-0 验证（`tleap/frame0_validation.txt`，`validate_frame0.py`）

- FAH 58/60 原子全连通、无孤立原子、羧基 C 价态 3 ✓；HSA 582 残基 ✓
- 远端羧基锚定：O1D→ARG483(NE) 2.78 Å、ARG346(NH2) 2.82 Å（两体系一致，= 晶体 MYR 的 2.80 Å 位姿）
- CPU 烟测：200 步最小化后 PE −1.084×10⁶ / −1.156×10⁶ kJ/mol（历史范围 −1.07×10⁶）

### run_md.py barostat bug（已修）

- `md/run_md.py` 原代码在 `Simulation()` **之后**才 `addForce(MonteCarloBarostat)`（注释却写 BEFORE）；OpenMM 8.5.2 实测不报错，但 barostat 不进入已建 Context（后加 Force 对已有 Context 无效）→ 原 NPT 加热/生产实为 NVT
- 修复：barostat 移至 `createSystem`/`add_restraints` 之后、`Simulation()` 之前；restart 分支也补上（原 restart 完全无 barostat）
- `common/scripts/run_md.py`（exp-A）此前已按 claude-Jun06 §3.2 修复，无需动

### MD 排队

- `md/launch_queue.sh` 已 setsid 启动（PID 2436125）：先 `while pgrep -f "common/scripts/run_md.py"` 等待 exp-A（GPU 2/3）结束，然后 GPU2=单酸 ×3、GPU3=二酸 ×3 串行，各 100 ns（50M 步）
- 日志：`md/queue.log` + `md/<system>/rep<i>/md_output.log`
- GPU 0/1（exp-D 保留）不使用

## 2026-07-20 — 6 replica MD 全部完成 + MM-GBSA/PBSA（单酸 vs 二酸）

### MD 完成
- **16:38** c18_monoacid rep3 完成（13.2 h）；**16:55** c18_diacid rep3 完成（13.2 h）
- queue.log："exp-C queue complete: 6 jobs finished"。2 体系 × 3 replica × 100 ns 全部就绪

### MM-GBSA/PBSA（沿用 exp-A 管线，`analysis/mmgbsa/`）

- 流程：cpptraj autoimage+strip（帧 502–1002 步长 4 → 126 帧/replica，即最后 50 ns）→ ante-MMPBSA 干拓扑（`-n ':583'`, mbondi2）→ MMPBSA.py 14.0 MPI 14 ranks/副本，GB(igb=5) + PB(inp=2 默认) + idecomp=1，单轨迹，**未含熵**
- rep1/rep2 四副本先行（~15:50 启动）；rep3 由 `watch_rep3.sh` 监视 MD 完成后自动抽取+启动（17:35 全部 exit 0）

**结果（3-replica mean±SD, kcal/mol）**：

| 体系 | ΔG_GB | ΔG_PB |
|---|---:|---:|
| 单酸 | −78.91 ± 2.17 | −37.15 ± 3.88 |
| 二酸 | −89.69 ± 4.02 | −55.64 ± 7.35 |
| **ΔΔG（二酸−单酸）** | **−10.78 ± 2.64** (Welch p=0.025) | **−18.49 ± 4.80** (p=0.030) |

**结论：二酸显著优于单酸，GB/PB 同向且 replica 间一致（3/3），设计逻辑成立。**

机制（分解 + 成对作用能 + 距离三重证据）：
- 远端羧基双盐桥锚定 ARG346+ARG483（2.72–2.75 Å，全部 6 轨迹稳定；裸静电 单酸 −249.8 / 二酸 −228.4）
- **二酸近端羧基与 ARG408 形成第二个稳定盐桥**（rep1/2/3 = 3.73/3.04/2.75 Å；裸静电 −97.8；分解 ARG408 贡献 单酸 ≈0 → 二酸 −7.9 GB / −8.5 PB）——这是 ΔΔG 的主要来源
- vdW 两体系完全相同；ΔΔG 全部来自"第二羧基静电增益（ΔEEL −99.6）− 去溶剂化代价（ΔE_GB +89.8）"的净正余额
- **勘正历史观察**：旧文档"近端羧基游离、距最近 Arg ~13 Å"在游离脂肪酸 @ FA3 体系不成立（3 replica 均锚定 ARG408）；该描述应只适用于脂肪酸接肽/连接子的构建

- 统计说明：副本内 SEM 经 n_eff 修正（0.6–2.2）；副本间 n=3，Welch t 为最低限度统计；二酸 rep3 最强（盐桥最紧）拉大 SD，幅度不宜当精确值
- 计算全在 CPU（峰值 ~84 ranks），未触碰 GPU；rep3 MM-GBSA 启动时 rep3 MD 已自然结束

---
*维护者：Kimi Code*
*最后更新：2026-07-20*

## 2026-07-20 — 链长系列构建与 MD 队列（U 形曲线验证）

### 对照轨迹确认

- c18_monoacid / c18_diacid ×3 replica ×100 ns 全部完成（7-19~7-20，GPU 2/3，exit 0，无 NaN；末帧 T≈310 K，PE≈−1.09×10⁶ kJ/mol，密度≈1.012）

### ⚠ 链长编号考据（重要，影响 U 形曲线横轴）

- 任务定义（化学正确）：Cn 二酸 = HOOC-(CH2)_{n-2}-COOH，C18 二酸 = 16 CH2（octadecanedioate, 54 atoms）
- **原 "c18" 对照实际为 18 CH2**：c18_diacid = 二十烷二酸（eicosanedioate，总碳 20，60 atoms）；c18_monoacid = 十九烷酸（nonadecanoate，总碳 19，58 atoms）。源自 2026-05-27 原构建脚本以"链碳数"而非"总碳数"命名
- 新系列按总碳数构建。c20_diacid（18 CH2, 60 atoms）与现有 c18_diacid 对照为**同一分子**（同协议，其 3 rep 轨迹可直接作为该链长数据点；本次 c20 ×3 rep 为同分子的额外 replica）
- **真正 C18 二酸（16 CH2, 54 atoms）尚无轨迹**——U 形曲线中心点缺失。补救：`python build_fa_fa3.py diacid 18 --force`（会覆盖 c18_diacid 相关文件，需先备份对照）或改名构建后补跑 3 rep
- 单酸/二酸对照的相对比较不受影响：两者链长相等（均 18 链碳），"二酸 vs 单酸"分离的正是额外羧基的效应

### 构建（build_fa_fa3.py v4 泛化：`diacid L [L ...]`，L=总碳数）

| 体系 | CH2 数 | FAH 原子 | 总原子 | Na+ | 净电荷 | O1D→ARG483(NE) | O1D→ARG346(NH2) |
|------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| c12_diacid | 10 | 36 | 86,404 | 16 | 0.0000 | 2.78 Å | 2.82 Å |
| c14_diacid | 12 | 42 | 86,410 | 16 | 0.0000 | 2.78 Å | 2.82 Å |
| c16_diacid | 14 | 48 | 86,416 | 16 | 0.0000 | 2.78 Å | 2.82 Å |
| c20_diacid | 18 | 60 | 86,428 | 16 | 0.0000 | 2.78 Å | 2.82 Å |
| c22_diacid | 20 | 66 | 86,434 | 16 | 0.0000 | 2.78 Å | 2.82 Å |

- 全部 tleap Errors=0；FAH 全连通无孤立原子；电荷精确 −2.000；远端羧基全部置于 MYR 1003 晶体位姿（O1D/O2D 与晶体 O 偏差 0.03/0.02 Å）；短链尾部自然埋于口袋浅处，未人为拉伸/折叠
- 验证记录：`tleap/frame0_validation_chains.txt`

### MD 队列（launch_queue_chains.sh，19:34 启动）

- GPU 0：c12 r1-3, c14 r1-2；GPU 1：c14 r3, c16 r1-3, c20 r1；GPU 2：c20 r2-3, c22 r1-3（各 5 任务串行，100 ns/个）；GPU 3 保留
- 启动 3 min 检查：首任务（c12r1/c14r3/c20r2）均通过 minimization（PE≈−1.345~−1.350×10⁶，与对照 −1.349×10⁶ 一致）进入 NPT 平衡，T≈310 K，PE≈−1.09×10⁶，密度≈1.01，~187 ns/d（每任务 ~13 h），无 NaN
- 预计 15 任务全部完成：~2.8 天后（2026-07-23 傍晚）

---
*维护者：Kimi Code*
*最后更新：2026-07-20*

## 2026-07-20(晚) — 补缺：c18true_diacid（真 C18 二酸）+ c16_monoacid（真 C16 单酸）

### 命名决策（确认前一条目的编号 bug）

- 现有 `c18_diacid`=C20（18×c3, 60 atoms）、`c18_monoacid`=C19（58 atoms），**文件一律不动**（MM-GBSA 分析已引用其 prmtop，覆盖会破坏可复现性）
- 新体系独立命名：`c18true_diacid`（HOOC-(CH2)16-COOH，真 C18 二酸 = 司美格鲁肽脂链）、`c16_monoacid`（CH3-(CH2)14-COOH，真 C16 单酸 = 利拉鲁肽链）
- 构建：`build_fa_fa3.py diacid 18 --label c18true_diacid`（脚本新增 --label，单体系时可用）；`build_fa_fa3.py mono 16`

### 拓扑与 frame-0 验证（tleap Errors=0，记录 `tleap/frame0_validation_wave2.txt`）

| 体系 | FAH 原子 | 总原子 | 电荷 | Na+ | O1D→ARG483(NE) | O1D→ARG346(NH2) | 最小化后 PE (CPU 200 步) |
|------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| c18true_diacid | 54 (22 重+32 H, 16×c3) | 86,422 | −2.000 | 16 | 2.78 Å | 2.82 Å | −1.233×10⁶ kJ/mol |
| c16_monoacid | 49 (18 重+31 H, 15 链碳) | 86,419 | −1.000 | 15 | 2.78 Å | 2.82 Å | −1.139×10⁶ kJ/mol |

- 两体系远端羧基均在 MYR 1003 晶体位姿（O 偏差 0.03/0.02 Å）、FAH 全连通无孤立原子、净电荷 0.0000、HSA 582 残基、无 clash（min FA-HSA 2.53 Å）

### 队列布局（不动正在运行的 launch_queue_chains.sh）

- **GPU 3**：`md/launch_queue_c18true.sh`（PID 1464359）——先等 exp-G pilot（`pgrep -f "exp-G"` 消失）；启动时 exp-G 无进程、GPU 3 空闲，故 19:56 立即开始 c18true_diacid rep1-3 串行。rep1 3 分钟检查：minimization OK（PE −1.344×10⁶）、NPT 平衡中 T≈307-311 K、PE≈−1.09×10⁶、186 ns/d、无 NaN。日志 `md/queue_c18true.log`
- **wave-2（GPU 0/1/2）**：`md/launch_queue_wave2.sh`（PID 1464356）——按 worker PID（1446459/1446460/1446461 → GPU0/1/2，映射取自 /proc/*/environ）分别等待 wave-1 对应卡结束后，在该卡跑 c16_monoacid rep1/2/3。等待条件校验 /proc/PID/cmdline 防 PID 复用。当前 3 个等待子进程均在等待态，wave-1 队列（PID 1446454）不受影响。日志 `md/queue_wave2.log`
- **注意**：若 exp-G 在 c18true 运行期间启动新任务，会与 GPU 3 竞争（exp-G 等待条件为启动时判定）

### 预计完成

- c18true_diacid ×3：~40 h → 2026-07-22 中午
- wave-1 15 任务：~2026-07-23 午后；c16_monoacid ×3 随后 ~13 h → 2026-07-24 凌晨
- 至此 U 形曲线数据点（真实总碳数）：C12, C14, C16, C18(=c18true), C20(=旧 c18_diacid + c20_diacid 共 6 rep), C22；单酸对照：C16(=c16_monoacid), C19(=旧 c18_monoacid)

---
*维护者：Kimi Code*
*最后更新：2026-07-20*

## 2026-07-23 — 链长系列 wave-1 MM-GBSA/PBSA 完成（c12–c22 二酸 × 3 rep）

wave-1 15 条轨迹（2026-07-23 午后全部跑完，无 NaN）→ 按 exp-C 既有协议算完 MM-GBSA（`analysis/mmgbsa/chain-series/`，详见 RESULTS-chain-series.md）。

**协议**：与 c18_monoacid/c18_diacid 完全相同——cpptraj autoimage+strip（帧 502–1002 步长 4 → 126 帧/rep = 最后 50 ns）、ante-MMPBSA `-n ':583'`（mbondi2）、MMPBSA.py 14.0 MPI 14 ranks/case、GB igb=5 saltcon=0.1 + PB 默认 inp=2 istrng=0.1、idecomp=1、单轨迹、无熵。15 case 分 3 批（rep1/2/3，batch_watch.sh 自动接续），全部 exit 0。

**结果（真实总碳数标注；3-rep mean±SD, kcal/mol）**：

| 体系 | C | ΔG_GB | ΔG_PB |
|---|---:|---:|---:|
| c12_diacid | 12 | −71.30 ± 9.52 | −57.01 ± 11.53 |
| c14_diacid | 14 | −80.87 ± 2.56 | −57.85 ± 4.77 |
| c16_diacid | 16 | −76.22 ± 2.57 | −45.83 ± 5.86 |
| c20_diacid | 20 | −88.32 ± 5.75 | −54.05 ± 5.40 |
| c22_diacid | 22 | −93.23 ± 1.03 | −51.98 ± 0.99 |

- **交叉验证**：新 c20_diacid vs 旧 "c18_diacid"（实为同一 C20 分子，FAH 原子名/电荷完全一致）旧值 GB −89.69±4.02 / PB −55.64±7.35 —— 差 ≤1.6 kcal/mol，管线自洽 ✓
- 形状（暂不画图）：GB 浅锯齿——c12 最弱、c16 局部低点、c20/c22 最强；PB 下 c12/c14 与 c20/c22 相当、c16 最弱。两方法一致点：c16 为局部低点
- 机制：远端双盐桥（ARG346/483）全链长稳定（2.72–2.76 Å）；近端 COO–ARG408 第二盐桥同样全链长存在（2.76–2.90 Å）；**FAH 贡献随链长单调增强（−13.4→−29.6 GB）**，是长链更强的主因
- **c12 rep3 锚点松散**（远端–ARG483/346 = 3.76/3.43 Å，配体未离袋）→ 该 rep 偏弱（GB −60.4），c12 SD 大即源于此；短链锚定弱本身即链长效应，如实保留
- 待办：c18true_diacid（真 C18）与 c16_monoacid（真 C16 单酸）MD 完成后按同协议补跑，再画 U 形曲线
- 全程 CPU（每批 70 ranks），GPU 任务未受影响

---
*维护者：Kimi Code*
*最后更新：2026-07-23*

## 2026-07-24 — 链长系列 wave-2 MM-GBSA 补缺 + U 形曲线（9 数据点）

c16_monoacid rep1/2/3（MD 全部完成）与 c18true_diacid rep1/2（rep3 仍在 MD，~43/100 ns @ 06:51）按与 wave-1 **逐字节相同**的 mmpbsa.in（GB igb=5 + PB inp=2 + decomp，126 帧/rep = 最后 50 ns，单轨迹，无熵）补跑完成：5 case × 14 ranks 一批（06:57 启动，08:17 全部 exit 0，watch_wave2.sh 监视）。产物汇总于 `analysis/mmgbsa/chain-series/`（RESULTS-chain-series.md、`U_curve.png`、`ucurve.py`、`ucurve_points.json`）。

**9 数据点汇总（真实总碳数；mean±SD, kcal/mol）**：

| 系列 | C | 体系 | n | ΔG_GB | ΔG_PB |
|---|---:|---|---:|---:|---:|
| 二酸 | 12 | c12_diacid | 3 | −71.30 ± 9.52 | −57.01 ± 11.53 |
| 二酸 | 14 | c14_diacid | 3 | −80.87 ± 2.56 | −57.85 ± 4.77 |
| 二酸 | 16 | c16_diacid | 3 | −76.22 ± 2.57 | −45.83 ± 5.86 |
| 二酸 | 18 | **c18true_diacid** | 2* | **−87.40 ± 3.21** | **−56.10 ± 3.19** |
| 二酸 | 20 | c20_diacid | 3 | −88.32 ± 5.75 | −54.05 ± 5.40 |
| 二酸 | 20 | c18_diacid（旧=C20） | 3 | −89.69 ± 4.02 | −55.64 ± 7.35 |
| 二酸 | 22 | c22_diacid | 3 | −93.23 ± 1.03 | −51.98 ± 0.99 |
| 单酸 | 16 | **c16_monoacid** | 3 | **−70.47 ± 1.51** | **−35.38 ± 4.97** |
| 单酸 | 19 | c18_monoacid（旧=C19） | 3 | −78.91 ± 2.17 | −37.15 ± 3.88 |

\* c18true 暂 n=2，rep3 MD 完成后补算更新。

**结论**：C18 **近最优但非唯一最优**——GB 下与 C20 同处平台期（差 ≤2，噪声内），但 C22（−93.2±1.0）仍更强（差 5.8，待 rep3 再判）；PB 下最优在 C14，C18 属顶部簇（C12/C18/C20 差异 ≤3）。两方法一致：c16 二酸是局部低点、C18 已进入顶部簇。单酸对照：C16 上二酸优于单酸（GB −76.2 vs −70.5，PB −45.8 vs −35.4），与"二酸优于单酸"主线结论一致。

- 分解：c18true 的 FAH 贡献 −21.4（GB TDC），落在链长单调趋势 c16(−20.3) < c18(−21.4) < c20(−25.5) < c22(−29.6) 上；c16_mono 的 ARG408 贡献 +0.01（无近端羧基，符合预期），远端锚 ARG346/483 −14.3/−15.9
- 锚点校验：新 5 轨迹远端双盐桥全部稳定（2.71–2.75 Å）；c18true 近端–ARG408 rep1 3.39±1.03 Å（rep2 2.78 Å）
- 提醒：FA3 结合能只是链长选择的一个因素；C18 vs C22 的 GB 差（~6）在 MM-GBSA 误差范围内，生物学结论应等 rep3 并结合其他证据
- 全程 CPU（70 ranks），未动 GPU（c18true rep3 的 MD 在 GPU 上继续，不受影响）

---
*维护者：Kimi Code*
*最后更新：2026-07-24*
