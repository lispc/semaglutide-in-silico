# exp-B: K34R 突变的结构补偿与酰化控制验证

## 概述

**对应决策**：决策 2（Lys34→Arg34 实现定点酰化控制 + 意外获得额外受体互作）

**核心假说**：
1. **工艺必要性**：Lys34 若不突变为 Arg，会在 Lys26 脂化时也被非特异性脂化，导致肽-受体界面破坏
2. **结构红利**：Arg34 的胍基通过水介导氢键与 GLP-1R Glu27 形成额外静电互作（4ZGM 晶体结构揭示），增强受体结合

**实验类型**：常规分子动力学 (cMD) + 氢键网络分析

## 背景

司美格鲁肽的前体肽在 Lys26 处引入脂链以实现 HSA 结合。然而 GLP-1(7-37) 天然含有两个 Lys（Lys26 和 Lys34）。如果不加控制，脂化反应会同时发生在两个位点：

- **Lys26 脂化**：正确位点。脂链从 α-螺旋表面伸出，不影响 GLP-1R 结合
- **Lys34 脂化**：错误位点。Frimann 2023 的 MD 模拟显示，Lys34 处的 FA 链指向 ECD 外侧而非疏水 patch，导致结合力显著下降

Novo Nordisk 的解决方案是将 Lys34 突变为 Arg34：
- 保持位置 34 的正电荷（K→R 保守替换），维持与周围残基的静电环境
- 阻止非特异性脂化（Arg 的胍基不被 NHS 活化酯酰化）
- 意外发现：Arg34 通过水分子与 GLP-1R Glu27 形成额外氢键网络（Lau 2015, Figure 8）

**关键实验证据**（Lau 2015）：
- 4ZGM 晶体结构（2.70 Å）明确显示 Arg34 胍基–H₂O–Glu27 氢键网络
- 类似物 6（无 K34R 突变）的 GLP-1R potency 为 7.1 pM，司美格鲁肽（含 K34R）为 6.2 pM——K34R 贡献 ~15% 的额外亲和力
- 类似物 3（K34 酰化）的 potency 显著下降

## 实验设计

### 三个体系

| 体系 | 肽序列（残基 7–37） | 脂链位置 | 目的 |
|------|---------------------|----------|------|
| B1 (K34) | H[Aib]EGTFTSDVSSYLEGQAAKEFIAWLV**K**GRG | 无 | 基线：K34 版的受体结合 |
| B2 (R34) | H[Aib]EGTFTSDVSSYLEGQAAKEFIAWLV**R**GRG | 无 | 核心对比：R34 的额外氢键 |
| B3 (K34-acyl) | H[Aib]EGTFTSDVSSYLEGQAAKEFIAWLV**K(C16)**GRG | Lys34-C16 | 验证：错误酰化的后果 |

所有体系均含 Aib8（从 exp-A 继承），受体为 GLP-1R ECD。

**为什么 B1/B2 无脂链**：K34R 的核心作用是控制酰化位点和提供额外受体互作。去除脂链（Lys26 也不脂化）可以**隔离 K34R 的纯结构效应**，不与脂链-HSA 效应混淆。脂链的影响将在 exp-C（HSA 结合）和 exp-D（Linker）中系统性研究。

### 受体与模板

| 组分 | 来源 | 说明 |
|------|------|------|
| GLP-1R ECD | PDB 3IOL / 4ZGM | 人 GLP-1R 胞外域（残基 24–145），含 α/β/α 三明治折叠 |
| 模板复合物 | 3IOL (2.10 Å) | GLP-1R ECD + 天然 GLP-1(7-37)，无脂链，WT 序列 |
| 对照模板 | 4ZGM (2.70 Å) | GLP-1R ECD + 司美格鲁肽（Aib8, R34, K26-酰化） |

**策略**：从 3IOL 出发（天然 WT 复合物），在肽上引入 Aib8（已有参数）+ K34R 突变。4ZGM 作为 Arg34–Glu27 氢键的参考坐标。

### 结构构建

1. 提取 3IOL 的 GLP-1R ECD + GLP-1(7-37) 坐标
2. 在肽位置 8 引入 Aib8（使用 exp-A 的 Aib 参数和 ParmEd 修改方案）
3. B1：保持 Lys34（无需修改）
4. B2：Lys34→Arg34 突变（tleap 或 ParmEd 直接改残基名+原子）
5. B3：在 B1 基础上，Lys34 NZ 连接 C16 脂肪酸（通过酰胺键，模拟 NHS 活化酯脂化产物）
6. 所有体系用 tleap 构建拓扑（ff14SB + GAFF2 + TIP3P），solvate + neutralize

### 对接/放置策略

3IOL 已有正确的肽-ECD 结合姿势，只需：
1. 将 Aib8 肽的 backbone 对齐到 3IOL 的 GLP-1 骨架
2. 局部能量最小化（肽 sidechain + ECD 界面残基）
3. 不从头对接——晶体结构提供了可靠的初始结合模式

### MD 方案

| 参数 | 设定 |
|------|------|
| 引擎 | OpenMM 8.5.1 |
| 力场 | ff14SB + GAFF2 + TIP3P |
| 非标准残基 | Aib8 (已有参数), Lys-C16 (需参数化) |
| 体系大小 | ~80–100k atoms（ECD ~14 kDa + 31-mer 肽 + 水） |
| 温度 | 310 K (Langevin) |
| 压力 | NVT production（同 exp-A） |
| 步长 | 2 fs (SHAKE on H-bonds) |
| 非键截断 | 10 Å, PME |

**第一阶段（验证）**：
- 3 systems × 1 replica × **200 ns**
- 验证体系稳定性、R34-Glu27 氢键形成

**第二阶段（生产）**：
- 3 systems × 3 replica × 300–500 ns
- 统计显著性和 replica 间变异

### 约束策略

相比 exp-A（DPP-4 催化腔需要约束防止长肽漂移），exp-B 的 GLP-1R ECD 结合是典型的 **α-螺旋–ECD 表面互作**：

1. **ECD backbone CA 弱约束**（10 kJ/mol/nm²）：防止受体整体平移
2. **肽 C 端 α-螺旋弱锚定**（残基 22–37 backbone, 50 kJ/mol/nm²）：维持螺旋在 ECD 疏水沟中的大致位置
3. **肽 N 端 + 中央区域不约束**（残基 7–21）：让 Aib8 区域和中央 loop 自由采样
4. **K34/R34 侧链 + Glu27 侧链完全自由**：这是实验的核心观察区

### 分析管线

| 分析 | 方法 | 预期 |
|------|------|------|
| 体系稳定性 | 骨架 RMSD（ECD Cα, 肽 backbone） | 200 ns 内收敛至 <3 Å |
| **R34-Glu27 氢键网络** | Arg34(NH1/NH2/NE) ↔ H₂O ↔ Glu27(OE1/OE2) 距离和 occupancy | B2 中 >50% occupancy, B1 中无 |
| **水分子桥接分析** | Arg34–Glu27 间桥接水分子数分布 | B2: 1–2 个桥接水稳定存在 |
| ECD 结合界面 RMSD | 肽残基 22–37 与 ECD 的界面 RMSD | B3 显著高于 B1/B2 |
| 肽 C 端螺旋稳定性 | 残基 22–37 的 helicity (DSSP) | B3 中螺旋可能部分解旋 |
| 接触频率 | Phe28/Leu32/Trp31 与 ECD 疏水 patch 的接触 | B3 中接触减少 |
| K34-C16 脂链取向 | 脂链末端与 ECD 的距离和角度 | B3: 脂链指向 ECD 外侧（非疏水沟） |

### 成功判据

- B2 中 Arg34–H₂O–Glu27 氢键网络 occupancy > 50%（相对于 B1 的 <10%）
- B3 中 ECD 结合界面 RMSD 显著高于 B1/B2（>1 Å 差异）
- B3 中 Phe28/Leu32 与 ECD 疏水 patch 接触数显著减少
- 体系无 catastrophic clash，PE < 10⁶ kJ/mol

## 软件与环境

| 工具 | 路径/环境 | 用途 |
|------|----------|------|
| Python | `gmx` env | 脚本运行 |
| OpenMM 8.5.1 | `gmx` env | 生产 MD |
| AmberTools 24.8 | `cgas-md` env | tleap, antechamber, parmchk2 |
| PyMOL | miniforge3 | 可视化 |
| mdtraj | `gmx` env | 轨迹分析 |

## 目录结构

```
exps/exp-B/
├── README.md          # 本文件
├── tasks.md           # 任务清单
├── exp-log.md         # 实验日志
├── structures/        # 输入的 PDB/结构文件
├── docking/           # 对接/结构准备脚本
├── tleap/             # tleap 拓扑构建
├── md/                # MD 输入/输出
│   ├── b1_k34/        # B1: Aib8, Lys34, 无脂链
│   ├── b2_r34/        # B2: Aib8, Arg34, 无脂链
│   └── b3_k34acyl/    # B3: Aib8, Lys34-C16
├── analysis/          # 分析脚本与结果
└── archive/           # 旧版本归档
```

## 参考

- Lau et al. (2015) *J. Med. Chem.* 58, 7370–7380 — Figure 8 (Arg34–Glu27 H-bond network), Table 1 (potency data)
- 4ZGM (RCSB) — GLP-1R ECD + semaglutide, 2.70 Å
- 3IOL (RCSB) — GLP-1R ECD + WT GLP-1, 2.10 Å
- Frimann et al. (2023) *J. Biomol. Struct. Dyn.* 41(11), 5007–5021 — Lys34 酰化 MD
- Knudsen & Lau (2019) *Front. Endocrinol.* 10, 155 — semaglutide discovery review
- best-practice.md — 计算实验规范

---

*创建日期: 2026-05-27*
*状态: 设计阶段，待确认*
