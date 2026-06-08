# Minimal Model 已知问题记录

## 问题 1: ECD C-terminal 断裂

### 现象
tleap 构建时报警:
```
There is a bond of 14.991 angstroms between C and N atoms:
.R<SER 129>.A<C 10> and .R<PRO 130>.A<N 1>
```

### 成因
- 7OR0 (cryo-EM 全受体结构) 中 residues 130–136 密度缺失，PDB 文件中不存在这些残基
- GLP-1R ECD 的真实边界是 **residues 1–128**，最后一个二硫键为 Cys85–Cys126
- Residues 129–145 是连接 ECD 和 7TM 的 **linker/stalk（铰链区）**，天然柔性大，cryo-EM 中常无法解析
- 当前脚本截断到 residue 145，把未解析的柔性 linker 也包了进来

### 文献依据
- PhD thesis (Frederiksen): "The N-terminal domain (amino acids 1 – 128)... stabilized by three disulfide bonds (Cys46–Cys71, Cys62–C104, Cys85–Cys126)"
- 7OR0 论文: ECD-focused refinement 只能解析到 linker 起始部分，密度在 linker 区域迅速衰减

### 可选方案
| 方案 | 操作 | 优点 | 缺点 |
|------|------|------|------|
| A. 截断到 residue 128 | 只保留 ECD 核心 | 避免引入未解析 loop；结构紧凑；与文献一致 | 丢失 linker，但 linker 不直接参与 peptide 结合 |
| B. Modeller 建模 loop | 补全 130–136 | 结构完整 | 需要额外计算；柔性 loop 可靠性低 |
| C. 保持现状 | 保留到 145，忽略断裂 | 最简单 | C-terminal 在 MD 中会异常摆动，可能扰动附近水层 |

### 选择
**方案 A**：将 ECD 截断点从 145 改为 **128**。这样最后一个完整残基是 CYS 126（参与二硫键），结构最稳定。

---

## 问题 2: Peptide–LNK 断开

### 现象
- complex_ecd.pdb 中 peptide C-terminal (GLY 126) 与 LNK 第一个原子距离 **19.01 Å**
- tleap 未自动形成 peptide-LNK 共价键
- MD 中 peptide 和 LNK 以独立分子运动

### 成因
- Semaglutide 的 LNK（γGlu-2×OEG-C18）通过 **Lys26 侧链** 连接到 peptide 骨架
- 7OR0 是 cryo-EM 全受体在去垢剂胶束中的结构，文献记载:
  > "only the linker Lys26-2xOEG moiety could be modeled... before the density merged into the lower-resolution regions of the detergent micelle"
- 7OR0 中的 LNK 坐标指向去垢剂胶束，而非与 peptide C 端连续；因此从 protein_final.pdb 提取时，peptide 和 LNK 之间出现 19 Å 间隙
- PDB 4ZGM（2.70 Å 晶体结构）有 semaglutide 骨架与 ECD 的高分辨复合物，但 **Lys26 是未酰化的**（无 LNK）

### 文献依据
- biorxiv (Zhang et al., 2021): 7OR0 中 linker 密度弱，只能建模 Lys26-2xOEG 部分，脂肪酸链伸入去垢剂胶束
- PDB 4ZGM 论文 (Lau et al., 2015): 晶体结构显示未酰化 semaglutide 骨架与 ECD 的结合模式

### 可选方案
| 方案 | 操作 | 优点 | 缺点 |
|------|------|------|------|
| A. 使用 4ZGM peptide + 手动连接 LNK | 用 4ZGM 高分辨坐标替换 7OR0 peptide，再平移 LNK 使其 N 端靠近 peptide C 端，tleap 加 bond | 4ZGM 晶体结构准确；peptide-ECD 结合构象可靠 | 需要手动对齐和加 bond |
| B. OpenMM 谐振约束 | 用 CustomBondForce 强制约束 peptide C 和 LNK N 距离 | 无需改拓扑 | 力常数需要调试；非真实化学键 |
| C. 保持断开 | 让两者自由运动 | 最简单 | 不符合 semaglutide 真实结构；LNK 可能无法正确到达 HSA |

### 选择
**方案 A**：
1. 用 `exps/exp-B/structures/pep_4zgm_clean.pdb` 替换 7OR0 的 peptide（4ZGM 是晶体结构，构象更准确）
2. 将 LNK（7OR0 residue 127）平移，使其 N16 原子靠近 peptide C-terminal GLY 的 C 原子（目标距离 1.3 Å，肽键标准距离）
3. 在 tleap 中手动添加 `bond` 命令连接两者

---

## 修复计划

1. 修改 `build_ecd_model.py`:
   - ECD 截断: 145 → 128
   - Peptide 来源: 7OR0 chain P → 4ZGM `pep_4zgm_clean.pdb`
   - LNK 平移: 将 LNK 移动到 peptide C 端附近，形成合理肽键几何

2. 修改 `build_ecd.in`:
   - 添加 `bond` 命令连接 peptide C-terminal 和 LNK N-terminal

3. 重新运行 tleap → 最小化 → 生产 MD

---

## 相关文件
- `exps/exp-B/structures/4ZGM.pdb` — 4ZGM 原始晶体结构
- `exps/exp-B/structures/pep_4zgm_clean.pdb` — 4ZGM 提取的 peptide（residues 10–37，R34）
- `exps/exp-B/structures/ecd_clean.pdb` — 4ZGM 提取的 ECD（residues 29–128）
- `exps/exp-F/build/lnk_noh_zero.mol2` — LNK 参数文件
