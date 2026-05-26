# exp-A: DPP-4 对 Aib8 突变体的结合阻断验证

## 概述

**对应决策**：决策 1（Ala8→Aib8 阻断 DPP-4 酶解）

**核心假说**：Aib8 的双甲基通过空间位阻显著降低 GLP-1 N 端肽与 DPP-4 活性位点的结合亲和力，从而保护肽键免受切割。

**实验类型**：常规分子动力学 (cMD) + MM-PBSA 结合自由能分析

## 背景

天然 GLP-1 在血液中的半衰期仅 1–2 分钟，主要原因是 DPP-4（dipeptidyl peptidase-4）在 Ala8–Glu9 之间特异性切割。DPP-4 的催化三联体（Ser630, Asp708, His740）隐藏在由 α/β 水解酶折叠域和 β 螺旋桨结构域构成的催化腔内。底物 N 端第二残基（P2 位）的侧链嵌入 S1 疏水口袋——对天然 GLP-1，Ala8 的单个甲基完美契合。

司美格鲁肽将 Ala8 替换为 Aib（α-氨基异丁酸），其 Cα 上连接两个甲基。双甲基引入的空间位阻使 Aib8 无法适配 S1 口袋，且极大限制了主链二面角的旋转自由度。

**关键实验证据**（Lau 2015）：
- 类似物 6（Ala8 版，其他与司美格鲁肽相同）的体内 PK 与利拉鲁肽相当（~12 h）
- 司美格鲁肽（Aib8 版）的体内 PK 为 46+ h
- 这说明仅靠增强的 HSA 结合不足以实现周给药——DPP-4 保护同样关键

## 实验设计

### 体系构建

| 组分 | 来源 | 说明 |
|------|------|------|
| DPP-4 受体 | PDB 1NU6, chain B | 人源 DPP-4, 2.10 Å |
| 天然底物肽 | GLP-1(7-37): HAEGTFTSDVSSYLEGQAAKEFIAWLVKGRG | 31 残基，N 端 free amine, C 端 carboxyl |
| Aib8 突变肽 | GLP-1(7-37) Aib8: H[Aib]EGTFTSDVSSYLEGQAAKEFIAWLVKGRG | 第 2 位 Ala→Aib |

**为什么用全长肽**：
- 提供完整的 GLP-1 结构上下文——C 端 α-螺旋的存在影响 N 端在溶液中的构象系综
- DPP-4 切割 GLP-1 时识别的不只是 N 端 2–3 个残基，全长肽的 C 端与 DPP-4 表面存在额外接触
- 后续实验 B–F 都在全长肽背景下进行，保持一致

**对抗长肽漂移的策略**（来自 CJC-1295 Trap 5/6 教训）：
- 肽残基 7–12（N 端区域）的 backbone 施加 moderate harmonic restraint（100 kJ/mol/nm²），锚定在活性位点
- 残基 13–37（C 端尾）完全不约束——允许自由运动，但不拖走 N 端
- DPP-4 非活性位点区域 CA 弱约束（10 kJ/mol/nm²）
- 催化三联体 + S1 口袋 + 肽键区域不约束——WT 和 Aib8 的行为差异在此处自由表现

### 对接策略

1. 从 1NU6 chain B 提取 Diprotin A（IPI 三肽）在活性位点的坐标作为模板
2. 使用 AF3 或 Boltz-2 预测 GLP-1(7-37) 单体结构，取 N 端 7–12 残基对齐到 Diprotin A 骨架
3. 将全长肽的 N 端（7–12）手动对齐到活性位点，C 端 α-螺旋（22–37）放置在 DPP-4 表面
4. 使用 Rosetta FlexPepDock 进行约束性对接（CoordinateConstraint on DPP-4 catalytic residues + peptide N-terminal backbone）
5. 选择催化攻击距离（Ser630 OG → Ala8/Aib8 C）最接近理想值（~2.8–3.2 Å）且整体 clash 最少的 pose 作为 MD 起始结构

### MD 方案

| 参数 | 设定 |
|------|------|
| 引擎 | **OpenMM 8.5.1**（env: `gmx`） |
| 力场 | **ff19SB** + TIP3P |
| 非标准残基 Aib | GAFF2 + AM1-BCC 电荷（通过 AmberTools antechamber 参数化） |
| 体系大小 | ~110–140k atoms（DPP-4 ~160 kDa + 31-mer 肽 + 水 + 0.1 M NaCl） |
| 温度 | **310 K**（生理温度，Langevin thermostat） |
| 压力 | 1 bar（MonteCarlo barostat, NPT） |
| 步长 | 2 fs（SHAKE on H-bonds） |
| 非键截断 | 10 Å（LJ switching 8→10 Å），PME 长程静电 |

**第一阶段（验证）**：
- 2 systems (WT, Aib8) × 1 replica × **200 ns**
- 目的：验证体系稳定性、催化几何维持、力场参数合理性

**第二阶段（生产）**：
- 2 systems × 3 replica × 300–500 ns
- 在第一阶段确认 setup 无误后启动，不急于开始

### 约束策略

参考 CJC-1295 Trap 5/6（29 残基长肽从口袋漂移），全长 31 残基肽必须采取约束策略：

1. **DPP-4 backbone CA restraint**：对非活性位点区域的 CA 施加弱 harmonic restraint（10 kJ/mol/nm²），防止受体整体漂移
2. **肽 N 端区域 moderate restraint**：对残基 7–12（His7–Phe12）的 backbone 重原子（N/CA/C）施加 harmonic restraint（100 kJ/mol/nm²），锚定 N 端在活性位点附近。这模拟了 DPP-4 对底物 N 端的天然识别和锚定
3. **肽 C 端区域不约束**：残基 13–37 完全自由，允许在溶剂中运动，但不拖走受约束的 N 端
4. **关键催化残基免约束**：Ala8/Aib8 的侧链、Ser630（OG）、Asp708、His740 完全不受约束——WT 和 Aib8 的差异在此区域自由表现

### 分析管线

| 分析 | 方法 | 预期 |
|------|------|------|
| 体系稳定性 | 骨架 RMSD（DPP-4 Cα, 肽 backbone） | 200 ns 内收敛至 <3 Å |
| 催化攻击距离 | Ser630 OG → Ala8/Aib8 羰基碳距离分布 | WT: 集中在 ~3 Å; Aib8: 更宽、均值更大 |
| S1 口袋适配 | Ala8/Aib8 侧链与 S1 口袋残基（Tyr547, Trp629, Tyr631 等）的距离 | Aib8 存在持续的空间冲突 |
| 结合自由能 | MM-PBSA/MM-GBSA（最后 100 ns） | ΔG_bind(WT) 显著更负 |
| 能量分解 | Per-residue decomposition | Aib8 在 S1 口袋内的范德华排斥贡献升高 |
| 灵活性对比 | 肽 backbone RMSF | Aib8 在 Ala8 位置的 RMSF 降低（刚性化），但整体结合不如 WT 紧密 |

### 成功判据

- WT 体系中 Ser630 OG 到 Ala8 羰基碳的平均距离 ≤ 3.5 Å（维持催化攻击姿势）
- Aib8 体系中该距离显著大于 WT，或分布更宽
- MM-PBSA ΔΔG = ΔG(Aib8) − ΔG(WT) > 0（Aib8 结合更弱）
- 无 catastrophic clash 或能量爆炸（PE < 10⁶ kJ/mol）

## 软件与环境

| 工具 | 路径/环境 | 用途 |
|------|----------|------|
| Python | `/home/scroll/miniforge3/envs/gmx/bin/python3` | 脚本运行 |
| OpenMM 8.5.1 | `import openmm` (env: `gmx`) | 生产 MD |
| AmberTools 24.8 | `/home/scroll/miniforge3/envs/cgas-md/bin/` | tleap, antechamber, parmchk2 |
| PyMOL | `/home/scroll/miniforge3/bin/pymol` | 可视化 |
| Rosetta 2026.15 | `/home/scroll/miniforge3/envs/rosetta/bin/` | FlexPepDock 对接 |
| Boltz-2 | `/home/scroll/miniforge3/envs/boltz/bin/boltz` | 备选：预测肽-DPP4 复合物结构 |

## 目录结构

```
exps/exp-A/
├── README.md          # 本文件
├── tasks.md           # 任务清单
├── exp-log.md         # 实验日志
├── structures/        # 输入的 PDB/结构文件
├── docking/           # Rosetta FlexPepDock 输出
├── md/                # MD 输入/输出
│   ├── wt/            # WT (Ala8) 体系
│   └── aib8/          # Aib8 体系
├── analysis/          # 分析脚本与结果
└── archive/           # 旧版本归档
```

## 参考

- Lau et al. (2015) *J. Med. Chem.* 58, 7370–7380 — 原始 SAR 论文
- 1NU6 (RCSB) — 人 DPP-4 晶体结构, 2.10 Å, Diprotin A 共晶
- Frimann et al. (2023) *J. Biomol. Struct. Dyn.* 41(11), 5007–5021 — 酰化 GLP-1 + GLP-1R MD
- best-practice.md — 计算实验规范
- CJC-1295/AGENTS.md — 已知陷阱（特别是 Trap 5: 长肽口袋漂移）

---

*创建日期: 2026-05-26*
*状态: 设计阶段，待确认*
