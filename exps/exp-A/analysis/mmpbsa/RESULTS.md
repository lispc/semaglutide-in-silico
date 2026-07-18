# exp-A MM-PBSA/GBSA 结果（WT vs Aib8, 最后 100 ns）

**计算日期**：2026-07-19 · **工具**：MMPBSA.py 14.0 (AmberTools, env `cgas-md`) · **执行**：MPI 16 ranks/体系（vendored mpi4py 4.1.2 + 系统 OpenMPI 4.1.2）

## 结论摘要

| 量 | GB (igb=5) | PB (inp=2) |
|---|---:|---:|
| ΔG_bind WT (kcal/mol) | +357.91 ± 17.38 | +522.33 ± 25.72 |
| ΔG_bind Aib8 (kcal/mol) | +350.44 ± 19.27 | +496.52 ± 22.79 |
| **ΔΔG = Aib8 − WT** | **−7.48** (SEM 2.29) | **−25.81** (SEM 4.99) |
| ΔΔG（clash 修正后） | −4.85 | −23.18 |
| 配对 t 检验 | t=−3.27, p=0.001, n_eff=135 | t=−5.17, p<0.001, n_eff=41 |

**方向与 roadmap 预期相反**：ΔΔG < 0，即在这套轨迹与单轨迹 MM-PBSA/GBSA 下，Aib8 的表观结合能**更负（更有利）**，而非更弱。绝对值巨大且为正（+350/+520 kcal/mol），无物理意义——被起始对接 pose 的 Lys696–Phe28  clash（+487 kcal/mol vdW）主导，见下文"假象分析"。**未含熵贡献，仅供相对比较。**

重要解读（详见"解读与警告"）：
- 假说相关的**局部**信号方向正确：突变位点本身（P2，肽第 2 残基）Aib8 比 Ala8 不利 +3.1 (GB) / +3.4 (PB) kcal/mol，主要来自 Aib 在 S1 口袋的极性去溶剂化惩罚（+11.6 vs +5.3）。
- 但全局 ΔΔG 被带电残基静电项（His7 N 端、Ser630、Asp708 等）淹没：Aib8 的 N 端被推出催化位姿（催化距离 5.25 vs 3.70 Å）后形成了静电上更"舒服"的非生产性接触（Ser630 静电贡献 −14.5 vs −2.3）。
- MM-PBSA 能量有利 ≠ 催化位姿可行。此前 5/5 几何指标（催化距离、接触数、S1 距离等）均显示 Aib8 被推离活性位点，与本结果并不矛盾。

## 方法

### 体系与轨迹
- WT：拓扑 `../../tleap/wt.prmtop`，轨迹 `../../md/wt/wt_traj.dcd`（2003 帧 × 100 ps = 200 ns）
- Aib8：拓扑 `../../md/aib8/aib8_modified.prmtop`，轨迹 `../../md/aib8/aib8_traj.dcd`
- 取样：帧 1000–2000、步长 5 → **201 帧**（对应 100.0–200.0 ns，与原几何分析"最后 100 ns"一致）
- 残基划分（两体系相同）：受体 DPP-4 = 残基 1–728（11662 原子），配体 GLP-1(7-37) 肽 = 残基 729–759（465 原子；730 = Ala8/Aib8），760–776 Na⁺，777+ 水
- 编号映射：prmtop 残基号 = PDB 号 − 38（Y547→509, W629→591, S630→592, Y631→593, D708→670, H740→702）

### 预处理（cpptraj V6.24.0）
1. `autoimage`（以 DPP-4 为锚）→ `strip :WAT,Na+,Cl-` → 写 NetCDF 干轨迹（`nobox`，12127 原子）
   - 输入脚本：`wt/extract_dry.cpptraj`、`aib8/extract_dry.cpptraj`；产物 `wt/wt_dry_last100ns.nc`、`aib8/aib8_dry_last100ns.nc`
   - 项目文档预警的"Rhombohedral vs Truncated octahedron"盒子不匹配**未出现**：prmtop 与 DCD 均为 Orthorhombic，cpptraj 直接读取无警告
2. `ante-MMPBSA.py` 生成干拓扑（`-s ':WAT,Na+,Cl-' -n ':729-759' --radii=mbondi2`；注意 ante-MMPBSA 不允许同时给受体+配体 mask，只给配体 mask，受体取补集）
   - 产物：各体系 `com_dry.prmtop`（12127）= `rec_dry.prmtop`（11662）+ `lig_dry.prmtop`（465），mbondi2 半径
3. 完整性校验：干轨迹催化距离 Ser630 OG→Ala8/Aib8 C = **3.70 ± 0.19 Å (WT) / 5.25 ± 0.29 Å (Aib8)**，方向与既有分析一致；受体-配体 COM 距离稳定（~34 Å, σ≈0.3），无 PBC 撕裂
   - 注：exp-log 2026-05-27 记为 5.0/6.0 Å；本次用现行 `full_analysis.py` 的原子选择直接复核原始 DCD 得 3.71 Å，与本结果一致，旧数值按当前脚本不可复现（推测旧版本选择或帧范围不同）

### MMPBSA.py 输入（`wt/mmpbsa.in`、`aib8/mmpbsa.in`）
- 单轨迹协议（complex/rec/lig 取自同一帧），201 帧全用
- GB：`igb=5, saltcon=0.100`（mbondi2 半径经 prmtop 传入）
- PB：默认 `inp=2`（sander ipb=2，已在 mdin 中确认），`istrng=0.100, radiopt=0, prbrad=1.4`
- 分解：`idecomp=1, dec_verbose=1`，打印残基 `509,591-593,670,702,729-759`（S1 口袋 + 催化三联体 + 全肽）
- **未算 nmode 熵**（太慢；结果仅供 WT vs Aib8 相对比较）
- 运行：`run_full.sh`（`mpirun -np 16 MMPBSA.py.MPI ...`，GB+PB+decomp 一次完成，WT ~2.0 h、Aib8 ~2.1 h，纯 CPU）
- 日志：`wt/mmpbsa_run.log`、`aib8/mmpbsa_run.log`；原始输出 `FINAL_MMPBSA.dat/.csv`、`FINAL_DECOMP.dat/.csv`

### 统计
- 逐帧 ΔG 取自 `FINAL_MMPBSA.csv`，用 `common/lib/stats.py` 算自相关时间、n_eff、SEM(std/√n_eff) 与配对 correlated t-test
- 脚本：`parse_results.py`（→ `parsed_results.json`）、`perframe_stats.py`（→ `perframe_stats.json`）

## 能量组分（kcal/mol，Differences 块）

| 组分 | WT (GB) | Aib8 (GB) | Δ | WT (PB) | Aib8 (PB) | Δ |
|---|---:|---:|---:|---:|---:|---:|
| VDWAALS | 302.26 | 297.43 | −4.83 | 302.26 | 297.43 | −4.83 |
| EEL | −107.12 | −198.92 | −91.79 | −107.12 | −198.92 | −91.79 |
| E_GB / E_PB | 189.38 | 278.38 | +89.00 | 224.06 | 296.00 | +71.93 |
| ESURF/ENPOLAR+EDISPER | −26.60 | −26.45 | +0.15 | 103.13 | 102.01 | −1.12 |
| ΔG_gas | 195.14 | 98.52 | −96.62 | 195.14 | 98.52 | −96.62 |
| ΔG_solv | 162.77 | 251.92 | +89.15 | 327.19 | 398.00 | +70.82 |
| **ΔG total** | **357.91** | **350.44** | **−7.48** | **522.33** | **496.52** | **−25.81** |

逐帧统计（n=201）：GB τ≈1–3 帧、n_eff≈63–179；PB τ≈5–13 帧、n_eff≈16–41。PB 的自相关明显更强，SEM(PB)=4.99 相应更大。

## 关键残基分解（TDC 总贡献，kcal/mol；GB / PB）

| 残基 | WT | Aib8 | Δ(Aib8−WT) | 说明 |
|---|---:|---:|---:|---|
| P2 位 (Ala8/Aib8, res 730) | −5.08 / −4.77 | −1.96 / −1.39 | **+3.12 / +3.38** | **方向符合假说**：Aib 极性去溶剂化惩罚 +11.6 vs +5.3，压过略优的 vdW/静电 |
| Y547 (res 509) | −0.11 / −0.13 | −5.81 / −4.99 | −5.70 / −4.86 | Aib8 新静电接触 (−8.8 vs +0.3) |
| W629 (res 591) | −1.36 / −0.90 | −1.17 / −0.94 | +0.20 / −0.04 | 无差异 |
| S630 (res 592) | −1.92 / +0.17 | −7.94 / −5.00 | −6.02 / −5.18 | Aib8 静电 −14.5 vs −2.3（非生产性接触） |
| Y631 (res 593) | −2.04 / −1.50 | −3.08 / −1.76 | −1.04 / −0.26 | 略有利 |
| D708 (res 670) | +0.78 / +9.99 | +0.28 / +0.71 | −0.50 / −9.28 | PB 下差异大（静电/溶剂化重排） |
| H740 (res 702) | −1.97 / −2.50 | −1.80 / −2.39 | +0.17 / +0.11 | 无差异 |
| His7 (res 729) | +3.06 / −9.44 | +4.70 / −3.57 | +1.64 / +5.87 | 自由 N 端，静电幅值巨大（±180），涨落主导 |
| Phe28 (res 750) | +241.2 / +242.3 | +241.4 / +242.7 | +0.20 / +0.38 | **clash 残基**，两体系相同 → ΔΔG 中抵消 |

## 假象分析：Lys696–Phe28 clash

- 预处理校验发现 ΔVDWAALS ≈ +300 kcal/mol（结合时 vdW 升**高**），异常。定位：DPP-4 Lys696（prmtop 658）CB/CG 与肽 Phe28（prmtop 750）芳环重原子间距 2.27–2.30 Å（正常 ≥3.3 Å），14–18 个重原子对 <2.8 Å，**贯穿整个最后 100 ns**，且**在最小化起始结构中已经存在**（+259 kcal/mol）——源于手动对接 pose（C 端螺旋放在 DPP-4 表面时 Phe28 穿入 Lys696 侧链），MD 无法解开（芳环拓扑锁死）。
- 成对能量定量（cpptraj `energy :658,750`，201 帧）：vdW = +487.1 ± 16.5 (WT) / +486.5 ± 17.0 (Aib8)，静电 −11.8 / −13.9 kcal/mol。**两体系几乎完全相同** → 在 ΔΔG 中相互抵消，修正量仅 ΔΔE = −2.62 kcal/mol。
- 因此**绝对 ΔG 无意义**（被 +487 的 clash 主导为巨大正值），但 ΔΔG 的假象污染有限（修正后 GB −4.85 / PB −23.18）。

## 解读与警告

1. **未含熵贡献**（无 nmode），数值仅供 WT vs Aib8 相对比较。
2. 单轨迹协议 + 单 replica；PB 逐帧自相关强（n_eff 低至 16），误差棒偏乐观。
3. ΔΔG 全局方向（Aib8 更有利）与 roadmap 预期（≥ +5 kcal/mol 强信号）**相反**，且 GB/PB 一致。其来源是带电残基静电项的重排（ΔEEL = −91.8 被 ΔE_GB/PB = +89/+72 大部抵消），属于 MM-PBSA 最不可靠的部分；叠加对接 pose 的 clash 假象，**本计算不能作为"Aib8 增强结合"的证据**，只能作为"单轨迹 MM-PBSA 在该体系上无法给出假说方向信号"的否定性记录。
4. 与既有几何结论的一致性解读：Aib8 的 N 端被挤出 S1 催化几何（距离 5.25 vs 3.70 Å、接触更少），能量上却在非生产性位姿找到静电补偿——DPP-4 抗性来自**无法占据催化位姿**，而非全局亲和力下降。建议后续用催化几何约束下的局部能量（如仅 N 端残基 729–740 的 decomp 子集和）或多 replica + 熵估算再检验。
5. 已知小问题不影响结果：aib8 prmtop 中 Aib 的 CB2 仍名 HA（仅命名）；`run_md.py` 打印显示肽 N 端 backbone 约束实际匹配 0 原子（`launch.log`），即生产中肽仅靠 DPP-4 CA 弱约束间接定位。

## 文件清单

```
analysis/mmpbsa/
├── run_full.sh              # 全量运行脚本（MPI）
├── parse_results.py         # FINAL_*.dat 解析 → parsed_results.json
├── perframe_stats.py        # 逐帧统计 → perframe_stats.json
├── vendor/                  # vendored mpi4py 4.1.2（不改动 conda env）
├── wt/  : mmpbsa.in, extract_dry.cpptraj, com/rec/lig_dry.prmtop,
│          wt_dry_last100ns.nc, FINAL_MMPBSA.dat/.csv, FINAL_DECOMP.dat/.csv,
│          mmpbsa_run.log, pair_658_750.dat (clash 对能量), test_*
└── aib8/: 同上（aib8_dry_last100ns.nc 等）
```
