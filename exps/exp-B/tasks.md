# exp-B 任务清单

> 状态: □ 未开始 | ● 进行中 | ✓ 完成 | ✗ 跳过

---

## Phase 0: 结构准备 (2–3 h)

- [ ] **0.0** 下载 3IOL 和 4ZGM PDB，提取 GLP-1R ECD + GLP-1
- [ ] **0.1** 构建三条肽的初始坐标
  - B1 (Aib8, Lys34)：从 3IOL GLP-1 出发，A8→Aib，保持 K34
  - B2 (Aib8, Arg34)：从 3IOL GLP-1 出发，A8→Aib，K34→R34
  - B3 (Aib8, K34-C16)：从 B1 出发，K34 NZ 接 C16 脂肪酸
- [ ] **0.2** 肽-受体复合物组装（对齐到 3IOL/4ZGM）
- [ ] **0.3** tleap 拓扑构建 (ff14SB + GAFF2 + TIP3P)
  - B1: WT 类似（只需 Aib8）
  - B2: K34R 突变（Amber 原生支持）
  - B3: Lys-C16 需参数化（antechamber + parmchk2）
- [ ] **0.4** Frame-0 验证（clash check, 关键距离, 初始 PE）
- [ ] **0.5** 能量最小化 + 短 NVT 平衡 (1 ns)

## Phase 1: 200 ns 验证 MD (1.5 天)

- [ ] **1.0** B1 (K34) 200 ns OpenMM，GPU 0
- [ ] **1.1** B2 (R34) 200 ns OpenMM，GPU 1
- [ ] **1.2** B3 (K34-acyl) 200 ns OpenMM，GPU 2
- [ ] **1.3** 监测前三小时（~50 ns）：RMSD 收敛、无 NaN、关键距离稳定

## Phase 2: 初步分析 (1 h)

- [ ] **2.0** 体系稳定性检查（ECD RMSD, 肽 RMSD）
- [ ] **2.1** Arg34–Glu27 距离时序（B1 vs B2）
- [ ] **2.2** 水分子桥接分析（B1 vs B2 的 R34/K34–Glu27 间水分子 occupancy）
- [ ] **2.3** B3 K34-C16 脂链取向 vs ECD 表面距离
- [ ] **2.4** 决策：R34 氢键网络是否形成？B3 脂链是否破坏界面？

## Phase 3: 全量分析 (1 h)

- [ ] **3.0** 取最后 100 ns 做全量结构分析
- [ ] **3.1** 氢键 occupancy 定量（Arg34–H₂O–Glu27）
- [ ] **3.2** ECD 结合界面 RMSD 时序对比（三体系）
- [ ] **3.3** Phe28/Leu32/Trp31 与 ECD 疏水 patch 接触数
- [ ] **3.4** 肽 C 端螺旋 helicity (DSSP) 对比
- [ ] **3.5** K34-C16 脂链 RMSF 和取向分布（B3）
- [ ] **3.6** 综合图表（6-panel: RMSD, R34-Glu27 dist, water bridge, contacts, helicity, lipid orientation）

## Phase 4: 决策门

- [ ] **4.0** 评估实验核心假说是否成立
- [ ] **4.1** 决定是否需要多 replica（Phase 5）
- [ ] **4.2** 记录结论到 exp-log.md

## Phase 5: 多 replica 生产 (2–3 天, 视需要)

- [ ] **5.0** B1 × 2 additional replica × 300 ns
- [ ] **5.1** B2 × 2 additional replica × 300 ns
- [ ] **5.2** B3 × 2 additional replica × 300 ns

## Phase 6: 最终分析 (2 h, 视需要)

- [ ] **6.0** Replica 间统计（correlated t-test, 报告 CV）
- [ ] **6.1** 综合报告（含所有 replica 的统计显著性）
- [ ] **6.2** 归档旧脚本和中间文件

---

## 已知风险与缓解

| 风险 | 缓解 |
|------|------|
| K34-C16 参数质量 | 用 antechamber + parmchk2 参数化 Nε-酰化 Lys，GAFF2 力场，参考 Liu 2025 的脂化方案 |
| 肽从 ECD 表面漂移 | C 端螺旋弱约束（50 kJ/mol/nm²），这是螺旋-ECD 疏水沟的标准互作——不像 exp-A 的 DPP-4 催化腔那样需要 N 端锚定 |
| tleap 拓扑构建失败 | 回退方案：ParmEd 直接修改 3IOL WT 拓扑（同 exp-A 策略） |
| B3 体系过大 | C16 链只增加 ~20 个重原子，体系大小增量 <5% |

---

*维护者：Claude Code*
*创建日期：2026-05-27*
