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
