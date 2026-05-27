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

### 待完成

- C18 monoacid 100 ns 验证完成后的初步分析
- 修复 diacid mol2 几何（O1P-O2P 距离问题）
- 构建其他脂链变体 (C16 monoacid, C12-C20 diacids)

---

*维护者：Claude Code*
*最后更新：2026-05-27*
