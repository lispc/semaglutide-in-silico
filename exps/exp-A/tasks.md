# exp-A 任务清单

> 标记：⬜ 待做 | 🔄 进行中 | ✅ 完成 | ❌ 阻塞/取消

---

## Phase 0: 结构准备

| ID | 优先级 | 状态 | 任务 | 说明 |
|----|-------|------|------|------|
| A-01 | P0 | ✅ | 下载并预处理 DPP-4 结构 (1NU6) | 提取 chain B，去除 NAG/溶剂，保留 Diprotin A 坐标作为参考 |
| A-02 | P0 | ✅ | 构建 GLP-1(7-37) WT 肽初始结构 | 31 残基全长，N 端 free amine, C 端 carboxyl。通过晶体结构 1NU6 + 序列比对手动构建 |
| A-03 | P0 | ✅ | 构建 GLP-1(7-37) Aib8 肽初始结构 | H[Aib]EGTFTSDVSSYLEGQAAKEFIAWLVKGRG，其他同 WT |
| A-04 | P0 | ⬜ | Aib 残基力场参数化（ff19SB 兼容） | 使用 antechamber + parmchk2 (GAFF2 + AM1-BCC)。非标准残基需与 ff19SB 蛋白力场兼容。输出到 `common/params/` |
| A-05 | P1 | ✅ | 验证 Aib 参数 | 最小化成功（PE 从 1.81×10¹⁵ → -1.73×10⁶ kJ/mol），200 ns MD 无 crash。实际用 ff14SB 非 ff19SB |
| A-06 | P1 | ⬜ | 验证 PDB 格式 | **阻塞**：`common/lib/pdb_validator.py` 仍未实现（roadmap Phase 0.6 承诺，至今未交付） |

## Phase 1: 对接

| ID | 优先级 | 状态 | 任务 | 说明 |
|----|-------|------|------|------|
| A-07 | P0 | ❌ | 准备 Rosetta FlexPepDock 输入 | **取消**：未使用 Rosetta 对接。肽通过 Kabsch 对齐 + 手动放置到活性位点 |
| A-08 | P0 | ❌ | 对接 WT 肽到 DPP-4 活性位点 | **取消**：同 A-07，未做 Rosetta 对接 |
| A-09 | P0 | ❌ | 对接 Aib8 肽到 DPP-4 活性位点 | **取消**：同 A-07，未做 Rosetta 对接 |
| A-10 | P1 | ❌ | 选择最佳 pose | **取消**：未做 Rosetta 对接，无 decoy 可选 |
| A-11 | P1 | ✅ | Frame-0 验证 | Clash check（无 VDW 重叠），最小化后 PE 正常（-1.73×10⁶ kJ/mol），催化距离 ~3 Å |

## Phase 2: MD 第一阶段（各 200 ns, 1 replica）

| ID | 优先级 | 状态 | 任务 | 说明 |
|----|-------|------|------|------|
| A-12 | P0 | ✅ | 构建 OpenMM MD 体系（WT） | 溶剂化 (TIP3P), 0.1 M NaCl, 中和电荷。体系 ~110–140k atoms |
| A-13 | P0 | ✅ | 构建 OpenMM MD 体系（Aib8） | 同上 |
| A-14 | P0 | ✅ | 编写通用 OpenMM MD 脚本 | `common/scripts/run_md.py` 已创建。ff14SB, 310 K, NVT production, DPP-4 CA 弱约束 + 肽 N-term BB 约束 |
| A-15 | P0 | ✅ | 运行 WT 200 ns（全长肽） | 200 ns OpenMM 完成（~178 ns/d），T=310K 稳定，催化距离 5.0±0.3 Å |
| A-16 | P0 | ✅ | 运行 Aib8 200 ns（全长肽） | 200 ns OpenMM 完成，催化距离 6.0±0.9 Å（比 WT 远 +1.0 Å） |
| A-17 | P1 | ✅ | 监控 MD 进程 | 22 ns / 47 ns / 200 ns 三次检查，全部记录到 exp-log |

## Phase 3: 第一阶段分析

| ID | 优先级 | 状态 | 任务 | 说明 |
|----|-------|------|------|------|
| A-18 | P0 | ✅ | 计算两个体系的骨架 RMSD 时序 | DPP-4 RMSD: WT 1.1±0.1 Å, Aib8 1.3±0.1 Å。肽 RMSD: WT 2.6±0.4 Å, Aib8 1.2±0.4 Å |
| A-19 | P0 | ✅ | 分析催化攻击距离分布 | Ser630 OG→Ala8/Aib8 C: WT 5.0±0.3 Å, Aib8 6.0±0.9 Å。**注意：当前只报了 mean±std，未做自相关校正/有效样本量（best-practice §32）** |
| A-20 | P0 | ✅ | 分析 S1 口袋空间冲突 | CB→Trp629: WT 11.3±0.4 Å, Aib8 13.0±0.4 Å（+1.7 Å，方向正确）。接触数: WT 1237 > Aib8 1188（-49） |
| A-21 | P1 | ⬜ | MM-PBSA/MM-GBSA 结合自由能 | **待做**：200 ns 数据就绪，但尚未计算。需 strip 水/离子后用 sander/MMPBSA.py |
| A-22 | P1 | ⬜ | Per-residue 能量分解 | **待做**：依赖 A-21 MM-PBSA |
| A-23 | P1 | ✅ | 肽 backbone RMSF 对比 | 已包含在 full_analysis.py 中 |
| A-23b | P1 | ⬜ | C 端漂移检查 | **待做**：未明确分析 |

## Phase 4: 决策门

| ID | 优先级 | 状态 | 任务 | 说明 |
|----|-------|------|------|------|
| A-24 | P0 | 🔄 | **Phase 1 评审** | 200 ns 结果已出（exp-log 2026-05-27），方向符合预期。等待用户正式确认是否进入多 replica 生产 |
| A-25 | P0 | 🔄 | **决定：是否进入多 replica 生产** | 同 A-24，等待用户决策 |

## Phase 5: 生产（多 replica + 长模拟）

| ID | 优先级 | 状态 | 任务 | 说明 |
|----|-------|------|------|------|
| A-26 | P1 | ⬜ | WT × 2 additional replica × 300 ns | **阻塞**：等待 A-25 决策 |
| A-27 | P1 | ⬜ | Aib8 × 2 additional replica × 300 ns | **阻塞**：等待 A-25 决策 |
| A-28 | P1 | ⬜ | 整合所有 replica 的 MM-PBSA 分析 | **阻塞**：依赖 A-26/A-27 + A-21 |
| A-29 | P1 | ⬜ | （预留） | |

## Phase 6: 最终分析

| ID | 优先级 | 状态 | 任务 | 说明 |
|----|-------|------|------|------|
| A-30 | P1 | 🔄 | 生成最终图表 | 部分图表已生成（full_analysis.png），但缺少 MM-PBSA chart（依赖 A-21） |
| A-31 | P1 | ⬜ | 更新 README.md 结论 | **待做**：当前 README 仍写"Phase 1 评审后更新" |
| A-32 | P2 | ⬜ | 写入项目主 README/project_log | **待做** |

---

*P0 = 阻塞性，必须先完成*
*P1 = 正常优先级*
*P2 = 可延后*
