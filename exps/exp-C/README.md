# exp-C: HSA 脂链- linker 复合体结合系统性分析

## 概述

**对应决策**：决策 3（C18 二酸为最优脂链）+ 决策 4（γGlu-2×OEG linker 隔离脂链）

**核心假说**：
1. C18 二酸在 HSA FA3 位点提供最优结合自由能（链长 U 形曲线）
2. γGlu-2×OEG linker 不仅隔离脂链与 GLP-1R，更**主动定位脂链**使其近端羧基与 HSA 表面第二个碱性残基形成锚点
3. 二酸优于单酸的关键在于远端羧基的额外静电锚定

**实验类型**：常规分子动力学 (cMD) + MM-GBSA 结合自由能分析

## 背景

### 未升级前的发现 (2026-05-27~28)

游离脂肪酸（无 linker）的 MD 模拟揭示：
- **C18 单酸**：远端羧基稳定锚定在 ARG482（2.8 Å, 100 ns），烷基链自由摆动
- **C18 二酸**：远端羧基同样锚定在 ARG482，但**近端羧基游离于溶剂**（距最近 ARG 13 Å），无法形成双点锚定
- **结论**：单独的 C18 二酸不能实现 Liu 2025 描述的双羧基锚定模式——需要 linker 将近端羧基约束在 HSA 表面附近

### 升级动机

Liu 2025 模拟的是**完整司美格鲁肽**（肽 + γGlu-2×OEG + C18 二酸），ΔG_bind = −75.33 kcal/mol。我们的游离 FA 实验揭示了 linker 的核心作用：它不仅是被动的"隔离臂"，更是**主动的定位器**——将近端羧基精确放置在 HSA 表面第二个正电荷残基附近。因此，系统性比较必须在 linker 上下文中进行。

### 脂链探针（升级版）

所有脂链通过 γGlu-2×OEG linker 连接（即 linker-FA 复合体），模拟司美格鲁肽中 Lys26 侧链的化学环境：

| # | FA 类型 | 碳数 | 说明 |
|---|---------|------|------|
| 1 | C12 二酸 | 12 | 最短二酸 |
| 2 | C14 二酸 | 14 | MYR 等长 |
| 3 | C16 二酸 | 16 | 接近最优 |
| 4 | **C18 二酸** | 18 | **司美格鲁肽**（阳性对照） |
| 5 | C20 二酸 | 20 | 长链（预期 U 型弯折/熵罚） |
| 6 | C16 单酸 | 16 | 利拉鲁肽型（单酸对照） |
| 7 | C18 单酸 | 18 | 二酸 vs 单酸对比 |

每个 linker-FA 复合体的近端（通过酰胺键连接 linker）和远端（游离羧基）化学环境均与司美格鲁肽一致。

### Linker 结构

```
γGlu — OEG — OEG — CO — (CH2)n — COOH (二酸远端)
  |                    |
  NH2                 近端酰胺键
```

- **γGlu**（γ-谷氨酸）：HOOC-CH2-CH2-CH(NH2)-COOH，通过 γ-羧基与 OEG 形成酰胺键
- **OEG**（8-氨基-3,6-二氧杂辛酸）：-NH-CH2-CH2-O-CH2-CH2-CO-，重复 2 次
- 参数化：GAFF2 + 标准电荷

### HSA 结构准备

| 组分 | 来源 | 说明 |
|------|------|------|
| HSA | PDB 1E7G, chain A | 人血清白蛋白, 2.50 Å, 残基 3–584 |
| FA3 模板 | 1E7G MYR 1003 | FA3 位点，远端羧基距 ARG482 2.8 Å |
| Linker-FA | 自建 (GAFF2 mol2) | γGlu-2×OEG-Cn diacid/monoacid |

### MD 方案

| 参数 | 设定 |
|------|------|
| 引擎 | OpenMM 8.5.1 |
| 力场 | ff14SB (HSA) + GAFF2 (linker-FA) + TIP3P |
| 体系大小 | ~67k atoms |
| 温度 | 310 K (Langevin) |
| 步长 | 2 fs (SHAKE on H-bonds) |
| 非键截断 | 10 Å, PME |

**Phase 1（验证）**：C18 linker-FA ×3 replica × 100 ns
**Phase 2（生产）**：7 systems ×3 replica × 300 ns = 6.3 μs

### 分析管线

| 分析 | 方法 | 预期 |
|------|------|------|
| 体系稳定性 | HSA RMSD, linker-FA RMSD | <3 Å |
| 结合自由能 | MM-GBSA（最后 100 ns） | C18 二酸最负 |
| 双羧基锚定 | 远端 COO⁻→ARG482, 近端 COO⁻→第二 ARG/LYS | C18 二酸双点锚定成立 |
| Per-residue 能量分解 | MM-GBSA decomposition | linker 贡献 vs FA 贡献 |
| ΔG vs 链长曲线 | 7 体系对比 | U 形曲线，C18 最优 |
| 链长效应 | C20 末端 RMSF, C12 口袋占有率 | 长链弯折，短链脱离 |

### 成功判据

- C18 linker-二酸实现**双点静电锚定**（远端和近端羧基各与一个 ARG/LYS 形成盐桥）
- MM-GBSA ΔG 在 C18 处最优
- Linker 贡献的能量可忽略（仅起定位作用，非结合力来源）
- 与 Liu 2025 的 −75.33 kcal/mol 在方法学层面可比较

## 目录结构

```
exps/exp-C/
├── README.md          # 本文件
├── tasks.md           # 任务清单
├── exp-log.md         # 实验日志
├── structures/        # 输入结构文件
├── tleap/             # 拓扑构建
│   ├── linker_fa/     # linker-FA 复合体 mol2 文件
├── md/                # MD 输入/输出
│   ├── c18_monoacid/  # 游离 FA 结果（保留，基线参考）
│   ├── c18_diacid/    # 游离 FA 结果（保留）
│   └── linker_c18/    # linker-C18 二酸（新）
├── analysis/          # 分析脚本
└── archive/           # 归档
```

## 参考

- Liu et al. (2025) *J. Biomol. Struct. Dyn.* — cMD + MM-GBSA of semaglutide-HSA, FA3 primary site
- Lau et al. (2015) *J. Med. Chem.* 58, 7370–7380 — original SAR
- Knudsen & Lau (2019) *Front. Endocrinol.* 10, 155 — semaglutide discovery review
- Curry et al. (1998) *Nat. Struct. Biol.* 5, 827–835 — HSA FA binding sites
- 1E7G (RCSB) — HSA + myristic acid, 2.50 Å

---

*创建日期: 2026-05-27*
*最后更新: 2026-05-28 — 升级为 linker-脂链复合体*
