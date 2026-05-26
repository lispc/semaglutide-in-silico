# exp-A 任务清单

> 标记：⬜ 待做 | 🔄 进行中 | ✅ 完成 | ❌ 阻塞/取消

---

## Phase 0: 结构准备

| ID | 优先级 | 状态 | 任务 | 说明 |
|----|-------|------|------|------|
| A-01 | P0 | ⬜ | 下载并预处理 DPP-4 结构 (1NU6) | 提取 chain B，去除 NAG/溶剂，保留 Diprotin A 坐标作为参考 |
| A-02 | P0 | ⬜ | 构建 GLP-1(7-37) WT 肽初始结构 | 31 残基全长，N 端 free amine, C 端 carboxyl。使用 AF3 或 Boltz-2 预测单体结构 |
| A-03 | P0 | ⬜ | 构建 GLP-1(7-37) Aib8 肽初始结构 | H[Aib]EGTFTSDVSSYLEGQAAKEFIAWLVKGRG，其他同 WT |
| A-04 | P0 | ⬜ | Aib 残基力场参数化（ff19SB 兼容） | 使用 antechamber + parmchk2 (GAFF2 + AM1-BCC)。非标准残基需与 ff19SB 蛋白力场兼容。输出到 `common/params/` |
| A-05 | P1 | ⬜ | 验证 Aib 参数 | 真空最小化 + 短 MD (10 ps) 验证无 crash。ff19SB + GAFF2 混合体系测试 |
| A-06 | P1 | ⬜ | 验证 PDB 格式 | 运行 PDB validator（待从 common/lib/ 创建），检查列对齐（best-practice §6） |

## Phase 1: 对接

| ID | 优先级 | 状态 | 任务 | 说明 |
|----|-------|------|------|------|
| A-07 | P0 | ⬜ | 准备 Rosetta FlexPepDock 输入 | DPP-4 (chain B) + 肽 pre-packed |
| A-08 | P0 | ⬜ | 对接 WT 肽到 DPP-4 活性位点 | 以 Diprotin A 坐标指导初始位置。CoordinateConstraint on DPP-4 catalytic residues。10 decoys |
| A-09 | P0 | ⬜ | 对接 Aib8 肽到 DPP-4 活性位点 | 同上 |
| A-10 | P1 | ⬜ | 选择最佳 pose | 标准：Ser630 OG → Ala8/Aib8 C 距离 (~3 Å)，无 clash，低 Rosetta 能量 |
| A-11 | P1 | ⬜ | Frame-0 验证 | Clash check, 关键距离（Ser630–肽键, 催化三联体内部），初始能量（best-practice §35） |

## Phase 2: MD 第一阶段（各 200 ns, 1 replica）

| ID | 优先级 | 状态 | 任务 | 说明 |
|----|-------|------|------|------|
| A-12 | P0 | ⬜ | 构建 OpenMM MD 体系（WT） | 溶剂化 (TIP3P, 12 Å buffer), 0.1 M NaCl, 中和电荷 |
| A-13 | P0 | ⬜ | 构建 OpenMM MD 体系（Aib8） | 同上 |
| A-14 | P0 | ⬜ | 编写通用 OpenMM MD 脚本 | 放 `common/` 下。支持: ff19SB 力场, 310 K, NPT equilibration → NVT production, backbone restraint on 肽 7–12 + DPP-4 non-active-site, restr file 指定约束原子 |
| A-15 | P0 | ⬜ | 运行 WT 200 ns（全长肽） | 注意应对长肽漂移：确认 restraint 参数生效，检查加热/NPT 阶段催化距离无异常偏离 |
| A-16 | P0 | ⬜ | 运行 Aib8 200 ns（全长肽） | 同上 |
| A-17 | P1 | ⬜ | 监控 MD 进程 | 每 24h 检查 RMSD/催化距离/能量趋势，记录到 exp-log |

## Phase 3: 第一阶段分析

| ID | 优先级 | 状态 | 任务 | 说明 |
|----|-------|------|------|------|
| A-18 | P0 | ⬜ | 计算两个体系的骨架 RMSD 时序 | DPP-4 Cα + 肽 backbone（N 端 7–12 受约束区域 vs C 端 13–37 自由区域分开看） |
| A-19 | P0 | ⬜ | 分析催化攻击距离分布 | Ser630 OG → Ala8/Aib8 羰基碳。画直方图对比。这是最核心的指标 |
| A-20 | P0 | ⬜ | 分析 S1 口袋空间冲突 | Ala8/Aib8 侧链与 Tyr547, Trp629, Tyr631, Val656 等的距离分布。Aib8 额外的甲基预期产生持续 clash |
| A-21 | P1 | ⬜ | MM-PBSA/MM-GBSA 结合自由能 | 最后 100 ns，取 200 frames, 0.5 ns interval |
| A-22 | P1 | ⬜ | Per-residue 能量分解 | 关注 S1 口袋残基对 WT vs Aib8 的范德华贡献差异 |
| A-23 | P1 | ⬜ | 肽 backbone RMSF 对比 | N 端受约束区 vs C 端自由区；Ala8/Aib8 位置灵活性差异 |
| A-23b | P1 | ⬜ | C 端漂移检查 | 验证 C 端自由运动没有通过 allosteric 效应拖动 N 端 |

## Phase 4: 决策门

| ID | 优先级 | 状态 | 任务 | 说明 |
|----|-------|------|------|------|
| A-24 | P0 | ⬜ | **Phase 1 评审** | 与用户讨论 200 ns 结果。确认：(a) 催化几何维持，(b) WT vs Aib8 差异方向符合预期，(c) 力场参数无异常，(d) 体系稳定性达标 |
| A-25 | P0 | ⬜ | **决定：是否进入多 replica 生产** | 如果结果符合预期 → Phase 5。如果有问题 → 调试后重跑 Phase 2 |

## Phase 5: 生产（多 replica + 长模拟）

| ID | 优先级 | 状态 | 任务 | 说明 |
|----|-------|------|------|------|
| A-26 | P1 | ⬜ | WT × 2 additional replica × 300 ns | 共 3 replica, 每个 >300 ns |
| A-27 | P1 | ⬜ | Aib8 × 2 additional replica × 300 ns | 同上 |
| A-28 | P1 | ⬜ | 整合所有 replica 的 MM-PBSA 分析 | 报告均值 ± SD（replica 间） |
| A-29 | P1 | ⬜ | （预留） | |

## Phase 6: 最终分析

| ID | 优先级 | 状态 | 任务 | 说明 |
|----|-------|------|------|------|
| A-30 | P1 | ⬜ | 生成最终图表 | RMSD 时序、催化距离直方图、MM-PBSA bar chart、S1 口袋 3D 可视化 |
| A-31 | P1 | ⬜ | 更新 README.md 结论 | |
| A-32 | P2 | ⬜ | 写入项目主 README/project_log | |

---

*P0 = 阻塞性，必须先完成*
*P1 = 正常优先级*
*P2 = 可延后*
