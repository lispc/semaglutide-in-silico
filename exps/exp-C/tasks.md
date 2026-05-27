# exp-C 任务清单

> 状态: □ 未开始 | ● 进行中 | ✓ 完成 | ✗ 跳过

---

## Phase 0: 结构准备与脂链参数化

- [✓] **0.0** 下载 1E7G, 分析 FA 位点分布
- [✓] **0.1** 确定 FA3 位点对应的 MYR 分子编号 (MYR 1003, R348/R485 附近 2.8 Å)
- [✓] **0.2** 提取 FA3 位点 MYR 坐标作为脂链模板
- [●] **0.3** 构建 7 种脂链的 3D 结构
  - C18 monoacid ✓ (58 atoms with H, GAFF2, charge -1)
  - C16 monoacid: 待
  - C18 diacid: mol2 几何需修复（O1P-O2P 距离 1.2Å→2.2Å）
  - C16 diacid, C12 diacid, C14 diacid, C20 diacid: 待
- [✓] **0.4** 参数化（GAFF2 + 标准羧酸电荷，antechamber 不可用，手动构建 mol2）
- [✓] **0.5** C18 monoacid 定位至 FA3（以 MYR 1003 为模板，Kabsch 旋转+平移，0 clash）
- [✓] **0.6** tleap 构建 HSA + C18 monoacid 拓扑（86k atoms, 25.7k waters）
- [✓] **0.7** Frame-0 验证（min dist ~1.6 Å, 无严重 clash, 初始 PE 正常）

## Phase 1: C18 monoacid 3-replica 验证 MD（进行中）

- [●] **1.0** C18 monoacid ×3 replica × 100 ns: GPU 0/1/2, 242 ns/day
- [ ] **1.1** 初步分析（RMSD, FA3 结合稳定性, 关键距离）
- [ ] **1.2** 决策：方法学可靠？扩大到全脂链集合？

## Phase 2: 全脂链参数化与 MD

- [ ] **2.0** 修复 diacid mol2 几何（O1P-O2P 距离归一化至 2.2 Å）
- [ ] **2.1** 构建 C16 monoacid + 定位至 FA3
- [ ] **2.2** 构建 C12-C20 diacids + 定位至 FA3（从 C14 MYR 模板延伸/缩短碳链）
- [ ] **2.3** 批量 tleap 拓扑构建
- [ ] **2.4** C16 monoacid ×3 replica × 300 ns
- [ ] **2.5** C18 diacid ×3 replica × 300 ns（阳性对照）
- [ ] **2.6** C16 diacid ×3 replica × 300 ns
- [ ] **2.7** C12, C14, C20 diacids ×3 replica × 300 ns

## Phase 3: 全量分析

- [ ] **3.0** MM-GBSA 结合自由能（所有体系，最后 100 ns）
- [ ] **3.1** ΔG_bind vs 脂链长度曲线（核心图表）
- [ ] **3.2** Per-residue 能量分解（R348/R485 静电贡献）
- [ ] **3.3** 远端羧基静电贡献定量（diacid vs monoacid ΔΔG）
- [ ] **3.4** 脂链 RMSF（C20 弯折检测）
- [ ] **3.5** 口袋占有率分析
- [ ] **3.6** 综合图表

## 已知风险

| 风险 | 缓解 |
|------|------|
| antechamber 不可用 | 手动构建 mol2 + 标准 GAFF2 电荷，所有脂链使用一致方法保证可比性 |
| diacid O1P-O2P 几何 | 修复中：增加 perp 偏移使 O-O 距离从 1.2→2.2 Å |
| GAFF2 电荷精度 | 使用标准官能团电荷（COO-: C=+0.70/O=-0.80, CH2: C=-0.12/H=+0.06），一致性优先于精度 |
| 其他 MYR 分子参数缺失 | 已从 HSA PDB 中移除，仅保留 FA3 位点 |

---

*创建日期：2026-05-27*
*最后更新：2026-05-27 21:50*
