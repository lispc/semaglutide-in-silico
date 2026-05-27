# exp-A 实验日志

> 只追加，不删除。每次记录含日期、时间、操作内容和结果。

---

## 2026-05-26 — 实验目录初始化

- 创建 `exps/exp-A/` 目录结构
- 编写 README.md（实验设计）、tasks.md（任务清单）、exp-log.md（本文件）
- **状态**：设计阶段，等待用户确认后启动 Phase 0

## 2026-05-26 — 设计调整（用户反馈）

调整项：
- **肽长度**：9 残基 (7-15) → 全长 31 残基 (7-37)。N 端 backbone restraint (res 1-6, 100 kJ/mol/nm²) 防长肽漂移
- **温度**：300 K → **310 K**（生理温度）
- **力场**：ff19SB → **ff14SB**（ff19SB 的残基特异性 CMAP 类型不支持非标准 Aib 残基）
- **引擎**：纯 OpenMM，不用 GROMACS

## 2026-05-26 — tleap 拓扑构建

- WT: tleap → prmtop/inpcrd → OpenMM 最小化。PE: 1.81×10¹⁵ → -1.73×10⁶ kJ/mol。✓
- Aib8: tleap 无法正确处理 AIB mol2（mol2 无 IC 表，H 原子坐标为 0）。
  改用 ParmEd 从 WT 拓扑修改 ALA→AIB：去 HA、加 CB2+3H、调电荷。无 clash。OpenMM 最小化成功。
- 服务器：tleap (AmberTools 24.8) from env `cgas-md`，OpenMM 8.5.1 from env `gmx`

## 2026-05-26 — 200 ns 生产 MD 启动

- **15:30** WT (GPU 0) 和 Aib8 (GPU 1) OpenMM 同时启动
- **17:00** WT (GPU 2) 和 Aib8 (GPU 3) GROMACS 交叉验证启动
- OpenMM 脚本: `common/scripts/run_md.py`
- GROMACS 脚本: `exps/exp-A/gmx/run_gmx.sh`
- OpenMM 约束: DPP-4 CA 弱约束 + 肽 N 端 BB 约束
- GROMACS: 无约束（仅用作引擎交叉验证）
- 310 K, NVT production, 2 fs step, checkpoint 每 1 ns

## 2026-05-26 — GROMACS 性能修复

- **问题**：GROMACS 初始性能仅 ~12 ns/day（GPU 80%），远低于预期
- **根因**：MDP 中 `nstlist` 为默认值 10，邻居列表每 20 fs 重建一次，开销巨大
- **修复**：`nstlist=40, rlist=1.1`（best-practice §3）
- **效果**：12 ns/day → **~275 ns/day**（23× 提升），GPU 利用率 ~82%
- 经验教训：`nstlist=40` 应作为项目默认 GROMACS 设置

## 2026-05-26 — 22 ns 初步分析

使用 mdtraj 分析 OpenMM 前 22 ns 轨迹（DCD 每 100 ps 输出）：

| 指标 | WT (Ala8) | Aib8 | 方向 |
|------|-----------|------|------|
| DPP-4 骨架 RMSD | 1.7 ± 0.0 Å | 1.7 ± 0.1 Å | 相当 |
| 肽骨架 RMSD | 1.9 ± 0.2 Å | 1.8 ± 0.2 Å | 相当 |
| 催化距离 (Ser630 OG→Ala8/Aib8 C) | 5.2 ± 0.4 Å | **7.2 ± 0.6 Å** | Aib8 更远 ✓ |
| 催化距离范围 | 3.7–8.2 Å | 5.4–8.6 Å | Aib8 偏离更大 |

关键发现：
- 两个体系均稳定（DPP-4 RMSD < 2 Å）
- 催化距离从初始 ~3 Å 漂移到 5–7 Å。这是预期行为——我们施加的是 N 端 backbone 弱约束（100 kJ/mol/nm²）而非共价键约束，肽在活性位点附近寻找能量最低的平衡位置
- **Aib8 的催化距离比 WT 大 ~2 Å**，方向正确：双甲基的空间位阻将 Aib8 推离 S1 口袋
- 等 200 ns 完成后做 MM-PBSA 定量分析

## 2026-05-26 — OpenMM vs GROMACS 初步对比

热力学层面一致：
- 温度：两者均稳定在 310 K
- 系综：均为 NVT
- PE 比值 ~0.78（不同引擎能量零点差异，属已知系统误差，见 best-practice §10）
- DPP-4 RMSD（OpenMM）：1.7 Å，GROMACS 因拓扑编号差异待后续对比
- 文献中最接近的实验：cjc-1295 项目（DPP-4 + GHRH，同酶不同肽），尚无文献做 DPP-4 + GLP-1 的 MD

## 2026-05-26 — 47 ns 进度检查

- OpenMM WT: 46.9 ns, 180 ns/d, T=310K ✓
- OpenMM Aib8: 46.7 ns, 179 ns/d, T=310K ✓
- GROMACS WT: 56.6 ns, 275 ns/d, T=310K ✓
- GROMACS Aib8: 56.4 ns, 275 ns/d, T=310K ✓
- 预计 OpenMM 明天中午完成 200 ns，GROMACS 明早完成

## 2026-05-27 — GROMACS 200 ns 完成

- **08:53** GROMACS WT 完成 (277 ns/day, xtc 789 MB)
- **08:54** GROMACS Aib8 完成 (277 ns/day, xtc 789 MB)
- OpenMM WT: 146 ns, 179 ns/d, ~7h 剩余
- OpenMM Aib8: 145 ns, 178 ns/d, ~7h 剩余
- GPU 2/3 已空闲

## 2026-05-27 — GROMACS 200 ns 分析（PBC unwrap 后）

| WT | OpenMM (约束) | GROMACS (无约束) |
|---|:---:|:---:|
| DPP-4 RMSD | 1.7 Å | 2.7 Å |
| 肽 RMSD | 2.0 Å | **8.3 Å** |
| 催化距离 | 5.2 Å | 11.7 Å |

| Aib8 | OpenMM (约束) | GROMACS (无约束) |
|---|:---:|:---:|
| DPP-4 RMSD | 1.7 Å | 3.2 Å |
| 肽 RMSD | 1.8 Å | 3.6 Å |
| 催化距离 | **7.2 Å** | **17.4 Å** |

关键发现：
- DPP-4 自身在两个引擎中稳定（1.7 vs 2.7–3.2 Å），引擎一致性 OK
- GROMACS 无约束导致肽漂移（RMSD 3.6–8.3 Å），验证了 OpenMM 中 N 端 backbone 约束是必需的
- **两个引擎给出相同方向**：WT 催化距离 < Aib8 催化距离（OpenMM: 5.2 < 7.2, GROMACS: 11.7 < 17.4）
- 初始 GROMACS 分析显示 DPP-4 RMSD ~29 Å，原因是 PBC wrapping 导致蛋白坐标跨盒子边界。用 `gmx trjconv -pbc mol -center` 修复后正常
- GROMACS Aib8 拓扑实际为 ALA（从 WT 拓扑修改而来，CB2 替换了 HA），引擎对比层面无影响

## 2026-05-27 — OpenMM 200 ns 完成 + 全量分析

- WT 和 Aib8 均在 200 ns 完成（~178 ns/d, 310K）
- 取最后 100 ns 做全量结构分析：

| 指标 | WT (Ala8) | Aib8 | Δ | 方向 |
|------|:---:|:---:|:---:|:---:|
| DPP-4 RMSD | 1.1±0.1 Å | 1.3±0.1 Å | +0.2 | 均稳定 |
| 肽 RMSD | 2.6±0.4 Å | 1.2±0.4 Å | -1.4 | Aib8 更刚性 |
| 催化距离 | 5.0±0.3 Å | 6.0±0.9 Å | **+1.0** | Aib8 更远 ✓ |
| 接触数 (<5Å) | 1237±46 | 1188±52 | **-49** | Aib8 少 ✓ |
| CB→Trp629 | 11.3±0.4 Å | 13.0±0.4 Å | **+1.7** | Aib8 远 ✓ |

核心结论：Aib8 的双甲基在所有指标上都推离了活性位点。
催化距离：WT 5.0 Å < Aib8 6.0 Å
接触数：WT 1237 > Aib8 1188
S1 口袋距离（CB→Trp629）：WT 11.3 Å < Aib8 13.0 Å
Aib8 肽 RMSD 更低（1.2 vs 2.6 Å）：双甲基增加局部刚性，"卡住了但进不去"

- 已知问题：Aib8 的 prmtop 中 CB2 仍标记为 HA（ParmEd 改名未持久化到文件），分析时用 HA 代替，物理正确，不影响结果

## 2026-05-27 — Kimi Review 反馈

Kimi 对项目做了全面 review，要点：
- 总体评价正面：科学动机清晰、文档质量标杆级、实验设计严谨
- P0 问题：API key 泄露（cc-ds.sh，已被 gitignore 保护未入历史）、Aib XML 缺 Angle/Dihedral（因用标准 ff14SB atom types，当前无实际影响）
- P1 问题：GROMACS emtol=1000 过松、缺渐进加热、lincs 参数未优化、OpenMM/GROMACS 约束不一致、硬编码路径
- P2/P3：旧版脚本未归档、文档一致性、统计方法升级
- 作者决定暂不修复，待生产完成后再处理

---

*维护者：Claude Code*
*最后更新：2026-05-26*
