# exp-C: HSA 脂肪酸结合系统性分析

## 概述

**对应决策**：决策 3（从 C16 单酸到 C18 二酸的范式转移）

**核心假说**：C18 二酸在 HSA FA3 位点提供最优结合自由能。更短的二酸疏水埋藏不足，更长的二酸遭遇熵罚或 U 型弯折。远端羧基与 HSA 碱性残基（R348, R485 等）的静电锚定是二酸优于单酸的关键。

**实验类型**：常规分子动力学 (cMD) + MM-GBSA 结合自由能分析

## 背景

司美格鲁肽区别于利拉鲁肽的最关键创新是将 C16 单酸（利拉鲁肽）升级为 C18 二酸（司美格鲁肽）。Lau 2015 的系统 SAR（Tables 1-4）清晰展示了脂链长度和二酸/单酸选择对 HSA 亲和力的影响。

**已有工作**（Liu 2025, *J. Biomol. Struct. Dyn.*）：
- 用 1.5 μs cMD + Pep-GaMD + MM-GBSA 确定了 FA3 为司美格鲁肽-HSA 的首选结合位点
- ΔG_bind(FA3) = −75.33 kcal/mol，主要贡献来自静电和范德华
- 确定了关键残基：R348, R485, R410, K414 等

**我们没有做的**：系统性比较不同脂链长度/类型的结合自由能曲线。这是逻辑链中"C18 二酸是最优解"的定量证据。

## 实验设计

### 脂链探针

| # | 类型 | 链长 | 碳原子数 | 说明 |
|---|------|------|---------|------|
| 1 | C12 二酸 | C12 | 12 | 最短二酸 |
| 2 | C14 二酸 | C14 | 14 | 中间长度 |
| 3 | C16 二酸 | C16 | 16 | 接近最优 |
| 4 | **C18 二酸** | C18 | 18 | **司美格鲁肽**（阳性对照） |
| 5 | C20 二酸 | C20 | 20 | 长链（Lau 假说：U 型弯折 + 熵罚） |
| 6 | C16 单酸 | C16 | 16 | 利拉鲁肽型（单酸对照） |
| 7 | C18 单酸 | C18 | 18 | 二酸 vs 单酸对比 |

共 7 种脂链，直接比较二酸长度效应 + 远端羧基贡献。

### 体系构建

**简化策略**（Phase 1）：
- 从 1E7G（HSA + 8× myristic acid）出发
- 选择 FA3 位点的 MYR 分子作为模板
- 将 MYR 替换为不同的二酸/单酸链
- 脂链参数化：GAFF2 + AM1-BCC（同 Liu 2025 方案）

**完整策略**（Phase 2，视需要）：
- 将脂链通过 γGlu-2×OEG Linker 连接到 Aib8,Arg34-GLP-1 骨架
- 完整 semaglutide 类似物 + HSA 体系
- 体系更大（~80-90k atoms vs ~65k for HSA alone with FA）

### HSA 结构准备

| 组分 | 来源 | 说明 |
|------|------|------|
| HSA | PDB 1E7G, chain A | 人血清白蛋白, 2.50 Å, 残基 3–584 |
| FA3 模板 | 1E7G MYR 分子（靠近 FA3 位点） | 确定哪个 MYR 在 FA3 |
| 脂链配体 | 自建（antechamber + parmchk2） | GAFF2 + AM1-BCC 电荷 |

**FA 位点识别**：
- HSA FA1-FA7 结合位点由 Bhattacharya et al. (2000) 和 Curry et al. (1998) 确定
- FA3 位点：位于 HSA 亚域 IB-IIA 界面，关键残基包括 R348, R485, R410, K414
- Liu 2025 图 2 标注了 FA3 位点的确切位置

### 对接策略

- 从 Liu 2025 获取 FA3 位点中 C18 二酸的初始结合姿势
- 若无法获取，以 MYR 羧基和碳链位置为模板，延伸/缩短碳链
- 二酸的远端羧基朝向 R348/R485 区域（静电锚定位置）

### MD 方案

| 参数 | 设定 |
|------|------|
| 引擎 | OpenMM 8.5.1 |
| 力场 | ff14SB (HSA) + GAFF2 (脂链) + TIP3P |
| 体系大小 | ~65k atoms (HSA alone + FA + water) |
| 温度 | 310 K (Langevin) |
| 步长 | 2 fs (SHAKE on H-bonds) |
| 非键截断 | 10 Å, PME |

**第一阶段（验证）**：
- 7 systems × 1 replica（先只做 C16 二酸、C18 二酸、C16 单酸）× **100 ns**
- 目的：验证脂链参数、FA3 结合稳定性、方法可行性

**第二阶段（生产）**：
- 7 systems × 3 replica × 300 ns = 6.3 μs
- 全脂链长度曲线

### 约束策略

- HSA backbone CA 弱约束（10 kJ/mol/nm²）：防止受体漂移
- FA3 口袋以外区域约束
- 脂链完全自由：允许在 FA3 口袋内自由采样
- 远端羧基（二酸）完全自由

### 分析管线

| 分析 | 方法 | 预期 |
|------|------|------|
| 体系稳定性 | HSA RMSD, 脂链 RMSD | 200 ns 内收敛 <3 Å |
| 结合自由能 | MM-GBSA（最后 100 ns） | C18 二酸最负 |
| Per-residue 能量分解 | MM-GBSA decomposition | 远端羧基与 R348/R485 的贡献明确 |
| 脂链 RMSF | 碳原子 RMSF | C20/C22 末端碳 RMSF 升高（弯折） |
| 口袋占有率 | 脂链与 FA3 口袋残基距离 | 短链（C12）口袋占有率低 |
| 远端羧基距离 | 远端 COO⁻ → R348 NH1/NH2 距离 | C18 二酸最优距离 |

### 成功判据

- MM-GBSA ΔG_bind 在 C18 二酸处最优（最负），形成 U 形曲线
- C12 二酸 → C18 二酸：ΔG 逐步降低（疏水埋藏递增）
- C18 二酸 → C22 二酸：ΔG 不再降低或回升（熵罚/U 型弯折）
- 二酸 vs 单酸（C18）：ΔΔG ≥ 5 kcal/mol 偏向二酸
- 远端羧基对结合自由能的贡献明确（per-residue decomposition）

## 软件与环境

| 工具 | 路径/环境 | 用途 |
|------|----------|------|
| Python | `gmx` env | 脚本运行 |
| OpenMM 8.5.1 | `gmx` env | 生产 MD |
| AmberTools 24.8 | `cgas-md` env | antechamber, parmchk2, tleap |
| mdtraj | `gmx` env | 轨迹分析 |
| MMPBSA.py | `cgas-md` env | MM-GBSA 计算 |

## 目录结构

```
exps/exp-C/
├── README.md          # 本文件
├── tasks.md           # 任务清单
├── exp-log.md         # 实验日志
├── structures/        # 输入结构文件
├── tleap/             # 拓扑构建
├── md/                # MD 输入/输出
│   ├── c12_diacid/
│   ├── c14_diacid/
│   ├── c16_diacid/
│   ├── c18_diacid/    # 司美格鲁肽型
│   ├── c20_diacid/
│   ├── c16_monoacid/  # 利拉鲁肽型
│   └── c18_monoacid/
├── analysis/          # 分析脚本
└── archive/           # 归档
```

## 参考

- Liu et al. (2025) *J. Biomol. Struct. Dyn.* — cMD + MM-GBSA of semaglutide-HSA
- Lau et al. (2015) *J. Med. Chem.* 58, 7370–7380 — original SAR, Tables 1-4
- 1E7G (RCSB) — HSA + myristic acid, 2.50 Å
- Curry et al. (1998) *Nat. Struct. Biol.* 5, 827–835 — HSA FA binding sites
- Bhattacharya et al. (2000) *J. Biol. Chem.* 275, 38731–38738 — FA site mapping
- Knudsen & Lau (2019) *Front. Endocrinol.* 10, 155 — semaglutide discovery review

---

*创建日期: 2026-05-27*
*状态: 设计阶段*
