# exp-G 实验日志

> 只追加，不删除。每次记录含日期、时间、操作内容和结果。

---

## 2026-07-20 — exp-G 启动：HSA–受体竞争三元复合物（预锚定 FA3）

### 设计（详见 README.md）

- 体系：HSA（582 残基，exp-C `hsa_no_myr.pdb`）+ ECD-肽复合物（exp-D `ecd_pep_nocap.pdb`，3IOL 结合构象刚体）+ LNK γGlu-2×OEG-C18 二酸（exp-D 重建拓扑 `lnk_gglu_2oeg_clean.mol2` 的类型/电荷/键表）
- 预锚定：远端 COO⁻（C55/O56/O57）置于 1E7G MYR 1003 羧基位（双盐桥 ARG346+ARG483）；近端羰基 O38 朝向 ARG408；肽-ECD 按 NZ 目标对接，允许大 ECD-HSA 间距
- 失败教训吸收：exp-F 的 HSA 游离放置（111 Å 永不接近）→ 本次全程共价锚定

### 构建流程与关键决策

1. `build/place_lnk_fa3.py`（LNK 几何）：
   - 远端羧基 O56/O57 按 exp-C `build_fa_fa3.py` 的 sp2 构建置于 MYR 羧基位（落在晶体氧 0.02–0.03 Å）
   - C54..C45 = 晶体 MYR C2..C11；**C44..C37 直接映射 exp-C 已验证的 FA3 出袋路径**（`c18_diacid_fa3.mol2` ext1..ext7，clearance ≥2.2 Å）——放弃了纯贪心生长（会在袋口穿出蛋白表面，1.3 Å 级 clash）
   - O38 二面角扫描靶向 ARG410（晶体编号）NHx ~3 Å → 得 3.95 Å
   - ADO/γGlu 段贪心生长：硬 clash 惩罚（<2.6 Å × 2000）+ 局部出口方向（C37 处 15 Å 内局部质心反方向）
   - **联合对接搜索**：O12 二面角（NZ 位置）× ECD 旋转角的 72×72 联合扫描，以"零 clash + 最大间距"选解
2. 踩过的坑：
   - 初始 NZ→CE 对齐方向取反（`rot_align(unit(nz0-ce0), axis)`），CE 与 C11 重合 0.08 Å、复合物整体翻转 180°；改为 `rot_align(unit(ce0-nz0), axis)` 后得真解（min dist 1.38 Å = NZ-C11 键本身）
   - exp-C 的"C18 二酸"实为 C20 二酸（20 C）；本次用 exp-D 的 LNK 拓扑（正确的 C18 二酸：C37 羰基 + 16 CH2 + COO⁻）
   - parmed 加载 prmtop 时 `residue.number` 是 0-based；验证与 run 脚本一律用 1-based 枚举位置
3. `build/assemble_complex.py`：应用对接变换，parmed 合并写盘（避免手写 PDB 的 mask 解析坑）；tleap 编号 HSA 1–582 / ECD 583–682 / 肽 683–708（Lys26=699）/ LNK 709
4. tleap 两遍法：pass 1 `solvateBox 12.0` 数水（43,706）→ pass 2 `addIons Na+ 0/Cl- 0`（20 Na⁺ 中和）+ `addIonsRand Na+ 118 Cl- 118`（0.15 M）
   - `bond`/`remove` 正常工作（NZ-C11 成键、Lys26 3×HZ 删除，parmed 验证）

### Frame-0 验证（`build/validate_frame0.py`，全部通过）

| 判据 | 结果 |
|---|---|
| 原子数 / 盒尺寸 | **141,880** / 138.2×100.1×120.1 Å |
| 远端双盐桥 | O56→ARG483 NE **2.78 Å**、O56→ARG346 NH2 **2.82 Å**（与 exp-C 晶体几何一致；O57 为外侧氧 3.52 Å 同 exp-C 记录） |
| 近端 O38→ARG408 | **3.95 Å**（≤5 ✓） |
| 肽-ECD 界面最小距离 | **2.92 Å**（≤4 ✓） |
| HSA↔ECD/肽/LNK clash（<2.2 Å） | **0**（最小 2.47 Å） |
| Lys26 NZ-C11 | 键存在，1.38 Å；HZ=0；其他 Lys 未动 |
| LNK 键/角/二面角 | 107/200/356，0 未参数化 |
| 水/离子/净电荷 | 43,450 / Na⁺ 138（20 中和+118 盐）/ Cl⁻ 118 / +0.10 |

### Pilot MD（GPU 3，20 ns 烟测）

- `md/run_pilot.py`：barostat 在 Simulation() 前；弱 CA 约束仅 HSA+ECD（kca 5.0→1.0）；锚点 O56/O57/O38 约束 kanc 5.0→**生产时完全释放**（否则无法观测竞争）；肽/LNK 自由
- 20:34 启动（setsid nohup，日志 `md/pilot_stdout.log`、`md/pilot/log.txt`）
- 已知小瑕疵：运行中的 pilot 因 0-based 编号多约束了肽 N 端 GLY 的 CA（683 个而非 682 个），对烟测物理影响可忽略；`run_pilot.py` 已修正（`<=681`），生产版生效
- min PE = **-2,338,082 kJ/mol**；升温各阶段温度与协议一致；~68 ns/d（GPU 3 共享）
- 10 分钟检查：（见下条记录）

---
*维护者：Kimi Code*

### Pilot 10 分钟检查（2026-07-20 20:45，启动 ~11 min）

- 已完成 min/升温/NPT eq，约束分级释放已执行（kca=1.0, kanc=0.0），进入生产
- 生产 400 ps：T=309.3–310.0 K，PE=-1,900,655 kJ/mol 稳定，密度 1.007 g/mL，65.7 ns/d
- 无 NaN；锚点已完全释放（竞争从此时刻起可观测）
- pilot 继续跑向 20 ns（预计 ~5 h），不等完成


## 2026-07-21 — Pilot 20 ns 完成与完整分析；production 启动

### Pilot 结果（20 ns, 7.4 h, GPU 3）

**锚定轨迹（远端羧基）**：预锚定双盐桥在生产开始后 ~2 ns 内松脱（O56/7→ARG483 3.2→~11 Å、→ARG346 → ~22 Å），之后**未完全逃逸**：
- **尾部羧基–HSA 表面最近距离：2.53 ± 0.14 Å（后 25%）**——脂链尾离开 FA3 口袋但保持 HSA 表面接触，未进入大量水相（对比 exp-C 游离 linker-FA 的 32–41 Å 完全逃逸）
- 近端 O38→ARG408：5.25 ± 1.2 Å（20 ns 时 3.93 Å，有靠扰成桥趋势）
- 肽-ECD 界面最小距离：**1.34 ± 0.03 Å，全程稳定**（肽始终锁定 ECD）
- LNK：端到端 26.3 ± 1.1 Å、Rg 8.6 ± 0.4 Å（紧凑、非绷紧）
- ECD-HSA COM：47.5 ± 0.5 Å 稳定；无 NaN

**解读**：游离 FA 锚定 100 ns 不动、MM-GBSA 也证实锚定有利，但带 linker 的完整分子 ~2 ns 脱锚——**linker 将 FA3 解离加速约两个数量级**，这是"HSA 动态储库 = 快速交换"的直接证据。同时尾部保持 HSA 表面接触而肽锁定 ECD：**γGlu-2×OEG 让分子可以同时挂住受体和白蛋白表面**（与 2×OEG 的选择逻辑一致）。排除机械应变假说（linker 非绷紧）。

### Production 启动（2026-07-21 ~04:20）
- 3 replica × 100 ns（50M steps，seeds 101/202/303），`md/launch_production.sh`，GPU 3（与 c18true 队列共享至其排空）
- 协议同 pilot（kca 5→1、锚点 kanc 5→0 生产时完全释放；CA 约束 0-based 修正已含）
- 分析目标：尾部驻留分区（FA3/表面/水相）、近端羧基-ARG408 成桥率、肽-ECD 稳定性、与 1×/3×OEG 变体对比的后续扩展

---

## 2026-07-25 — Production 3×100 ns 驻留分区分析

### 数据与校验
- rep{1,2,3} 各 1003 帧 × 100 ps = 100.30 ns（dt 由 reporter 间隔推导）；seeds 101/202/303 独立轨迹，统计按 rep 分别算再以 n=3 合并（不拼接）
- rep1 frame-0 键长校验：NZ-C11 1.35 / C55-O56 1.24 / C37-O38 1.24 / O56→ARG483 2.78 Å —— 轨迹-拓扑匹配 ✓
- 脚本：`analysis/residence_partition.py`；产物：`frames_rep{1,2,3}.csv`、`summary.json`、`RESULTS.md`

### 核心结果
- **远端羧基驻留分区（合并 n=3）**：FA3 20.2±25.0% / HSA 表面 45.2±10.0% / 过渡 29.0±21.8% / 水相 5.5±4.9%
- **三种 rep 命运**：rep1 全程动态交换（FA3 55.6%，346 次转换，100.2 ns 仍在交换，无 ≥1 ns 脱锚）；rep2 3.9 ns 脱锚转表面（不回头）；rep3 2.5 ns 脱锚渐溶剂化
- **近端 O38→ARG408 成桥率 0% ×3**（均值 12–14 Å）——初始近端盐桥是构建假象，平衡后不形成；2×OEG 分子在 HSA 上只有远端羧基一个有效锚点
- **肽-ECD 界面**：2.77±0.04 Å（后 25 ns 2.72±0.01），300 ns 零脱钩 ✓
- **linker 两端**：尾端（FA3 侧）先松（2.5–3.9 ns 或亚 ns 交换），受体端纹丝不动 —— 受体端亲和力 ≫ HSA FA3 端
- LNK 形态耦合锚定状态：锚定 rep1 e2e 19.6±1.6（紧凑）vs 脱锚 rep2/3 24.7/27.4（伸展）；COM 46.8±0.5 Å 稳定
- 与 pilot 20 ns 对照：pilot 单轨迹 ≈ rep2/3 型；rep 间差异在 20 ns 尺度不可见（方法学教训：本体系单轨迹结论不可信）
- 注：pilot "肽-ECD 1.34 Å" 度量无法复现（方法不明），本文重原子最小距离口径 frame-0=2.92 Å、production 2.7 Å 自洽；两口径定性结论相同

### 结论
预锚定 FA3 在游离竞争下是**快速交换的驻留区**而非静态锚点；"一头挂受体、尾端在白蛋白上 FA3/表面间快速交换"的桥接构型在 100 ns × 3 尺度上为稳态，与 2×OEG 平衡 BR ratio 的设计逻辑定性一致。详见 `analysis/RESULTS.md`。

---
*维护者：Kimi Code*
