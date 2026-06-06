# 司美格鲁肽 in silico 项目 Review — Claude, 2026-06-06

> 本 review 基于对项目**全部 18 个 Markdown 文件**的完整阅读，并对其中的关键技术声明（barostat 时序、API key 状态、力场参数、代码状态）逐一回到源代码 / 文件系统核对。
>
> 覆盖文件：`roadmap.md`、`docs/best-practice.md`、两份文献文档（`司美格鲁肽药物发现模拟设计.md`、Lau 2015 全文）、exp-A/B/C/D 的 README + exp-log + tasks、以及三份历史 review（kimi-May26、kimi-Jun01、kimi-Jun01-response）。
>
> 与前两份 Kimi review 的区别：本次重点放在 **(1) 科学结论是否被数据支撑**、**(2) 文档声明与代码/数据实际状态是否一致**、**(3) 历史 review 承诺的修复是否真正落地**。Kimi 已充分覆盖工程债务清单，本文不重复罗列，只补充新发现或核实结论。

---

## 一、总体评价

这是一个**科学叙事极强、文档纪律罕见、自我纠错文化真实存在**的项目。把司美格鲁肽拆成"Aib8 → K34R → C18 二酸 → γGlu-2×OEG"四步决策链，并为每一步设计一个可证伪的 MD 实验，这种"决策链驱动"的结构在个人计算生物学项目里是顶级水准。`best-practice.md` 和 `exp-log.md` 的"只追加、诚实记录失败"原则在 exp-C 的负面结果（linker-FA 逃逸）中得到了真实检验——这一点必须高度肯定。

但本次 review 的核心结论是：**项目当前最大的风险已经从"工程债务"转移到"结论强度与文档可信度的落差"。** 具体表现为三层：

1. **exp-D 的最终结论（"定性吻合 Lau 2015 Table 3"）被数据严重高估** —— 变体间差异在误差棒以内、无统计检验、关键的结合自由能（MM-GBSA）根本没算，而 README 承诺的 500 ns×3 实际只跑了 ~100 ns。
2. **历史 review 中"已修复"的若干项实际并未落地** —— 最突出的是 `common/scripts/run_md.py` 的 barostat 时序 bug（Jun01-response 声称已修，代码里没改）。
3. **被反复点名的 P0 安全问题以新形态复发** —— 新出现的 `cc.sh`（2026-06-06，未被 gitignore）再次明文存放 token。

下面分三部分展开。

---

## 二、做得好的地方 ✅

### 2.1 科学叙事与文献定位（项目灵魂）

- **`roadmap.md` 的四步决策链 + A–F 实验映射**是整个项目的最大资产。每个实验都明确写出"为什么必须做""对应哪个决策""已有文献做了/没做什么"。`roadmap.md:53-58` 的文献定位表（Frimann 2023 / Liu 2025 / 两份博士论文）精准划出了学术空间，避免了重复造轮子。
- **两份文献文档质量高**：`司美格鲁肽药物发现模拟设计.md` 是一份带 SAR 数据表的高质量综述，可直接作为论文 Introduction 底稿；Lau 2015 全文在库，便于随时核对原始数据。

### 2.2 真实的自我纠错文化

- **exp-B 的暂停决策**（发现所有晶体/cryo-EM 结构中 ECD-肽分离 25–40 Å，FlexPepDock 无法桥接 → 果断暂停而非硬上）是 best-practice §29「知道什么时候停」的教科书示范。
- **exp-C 的假说证伪被完整记录**：原假说"linker 主动定位近端羧基实现双点锚定" → 实测"无论封端与否，linker-FA 都从 FA3 逃逸（32–41 Å）" → 给出物理解释（OEG 亲水性）。`exp-C/exp-log.md:171-186` 把四轮 NaN 诊断、ACE 封端证伪、最终结论一字不删地保留下来。这是真正的科学诚实。
- **Jun01-response 对 review 的逐条回应有理有据**：对"dead script"类指控（launch.sh 从未执行、build_complex.py 双键拓扑从未用于生产）的反驳是站得住脚的——我核对了 `launch.sh` 已删除，且生产用的 exp-C/exp-D `run_md.py` 中 barostat 确实在 context 首次使用前加入（见 §3.2）。

### 2.3 工程上的真实进步

- **`run_md.py` 从 exp-A 到 exp-C/D 的演进**：从硬编码 WT/Aib8 升级为 `--system/--replica/--gpu` 参数化，是明确的工程改善。
- **exp-C/exp-D 的 barostat 时序实际是对的**：尽管注释措辞含糊，`exps/exp-C/md/run_md.py:62,77` 和 `exps/exp-D/md/run_md.py:53,63` 中 `addForce(MonteCarloBarostat)` 都在 `context.setPositions()` + 首个 `step()` **之前**完成，因此 NPT 在这两个实验中是生效的。Kimi 对这两个文件的担忧不成立（问题只在 `common/scripts/run_md.py`，见 §3.2）。
- **exp-D 的 NZ-C 共价键修复链条值得称道**：从"tleap bond 静默失败 → NZ-C 距离 14 Å、变体差异被掩盖"这个**会直接毁掉结论**的问题，定位到 ParmEd `BondType(427.0, 1.38)` 修复，并用单 replica → 全量验证。这是一次高质量的 debug，且全程记录在 `exp-D/exp-log.md:99-167`。

---

## 三、做得不好甚至不对的地方 ❌

### 3.1 🔴 exp-D 的最终科学结论被数据严重高估（本次最重要的发现）

exp-D 的 README 成功判据和 commit message（`9de1e38 exp-D: final analysis — qualitative match to Lau 2015 Table 3`）都声称"计算复现了 Novo Nordisk 的 linker 设计逻辑"。但回到数据本身，这个结论**当前不成立**：

**(a) 变体间差异在误差棒以内，且无统计检验。** `exp-D/exp-log.md:204-210` 的最终结果：

| Variant | EC50 (pM) | CA RMSD (Å) | Tail-Prot (Å) |
|---------|-----------|-------------|---------------|
| No linker | 269 | 2.3±0.2 | 4.6±0.1 |
| gGlu | 9.9 | 2.2±0.2 | 4.5±0.7 |
| gGlu-1×OEG | 4.8 | 2.5±0.3 | 4.4±0.6 |
| gGlu-2×OEG | 6.2 | 2.5±0.2 | 3.8±0.2 |
| gGlu-3×OEG | 27.7 | 3.0±0.2 | 4.1±0.4 |

- CA RMSD 跨度 2.2–3.0 Å，相邻变体的 ±SD 大量重叠。把 3×OEG 的 3.0±0.2 解读为"过长 linker 的 entropic penalty"——但 1×OEG（EC50 最低=4.8）的 RMSD 是 2.5，反而高于 gGlu（EC50=9.9）的 2.2，趋势并不单调。
- Tail-Prot 的"最优"信号是 2×OEG 的 3.8 Å vs 其余 4.0–4.6 Å —— **0.2–0.8 Å 的差异，且这些 ±SD 是按 replica 还是按帧算的没有说明**，更没有 best-practice §11/§32 反复强调的 correlated t-test 或 replica 间 CV。按项目自己的铁律，这个差异**不能宣称显著**。

**(b) 实验数据本身就不支持"2×OEG 最优"这个待复现的结论。** Lau 2015 的 EC50（无 HSA）里，1×OEG=4.8 < 2×OEG=6.2 < gGlu=9.9。也就是说**纯受体效力上，2×OEG 并不是最优**——司美格鲁肽选 2×OEG 是 BR ratio（含白蛋白后的平衡）和综合性质的结果，而不是裸 EC50。用一个不在 GLP-1R 端单调的 EC50 序列去"定性吻合"几何指标，本身逻辑就有问题。`司美格鲁肽药物发现模拟设计.md:54` 已经正确指出 2×OEG 的优势在 BR ratio——但 exp-D 没有模拟 HSA 端，无法触及 BR ratio。

**(c) 承诺的核心定量指标（MM-GBSA）根本没算。** roadmap §实验D 和 best-practice 都把结合自由能作为核心产出，但 exp-D 全程只有几何指标（RMSD、距离），grep `exps/exp-D/` 无任何 MM-GBSA/MM-PBSA 计算。仅凭 0.x Å 的几何差异声称"复现设计逻辑"，证据强度远不够。

**(d) README 的方案与实际执行不符。** `exp-D/README.md:32-33` 写"NPT production (500ns)"、"3 replica × 500 ns = 7.5 μs"，但 `exp-D/exp-log.md:198` 实际是"13/15 轨迹可读 ~95–102 ns"，且 9/15 是从 checkpoint 续跑拼起来的。**实际模拟量是承诺的 ~1/5**，文档未在 README 中更正。

> **建议**：把 exp-D 的结论降级为"初步/不显著"。要支撑"复现 linker 设计逻辑"，至少需要：(1) 跑满 replica 并报告 replica 间 CV + correlated t-test；(2) 补 MM-GBSA；(3) 诚实说明裸 EC50 不单调、真正的最优性体现在 BR ratio（需要 HSA 端，即把 exp-C 和 exp-D 合并成完整司美格鲁肽体系）。在此之前，commit/README 里的"qualitative match"措辞应加上"within error / not statistically significant"的限定。

### 3.2 🔴 历史 review 承诺"已修"但实际未落地的项

**barostat 时序 bug（`common/scripts/run_md.py`）——声称已修，实际没改。**
- Jun01-response 表格明确写："Move `system.addForce(MonteCarloBarostat)` before `Simulation()` in both exp-C and exp-D `run_md.py`"，并把它列为"P0: Fix this week / Agreed"。
- 实际核对：`common/scripts/run_md.py:135` 是 `Simulation(...)`，而 `:169` 的 `addForce(MonteCarloBarostat)` 在其后、且在 minimization + 第一段 NVT heating 已经跑过之后才加入；全仓 `grep reinitializeContext` **无任何结果**。
- 也就是说 **exp-A 用的那个脚本的 NPT 加热段实际跑在 NVT 下**。所幸 exp-A 用 NVT production，影响仅限平衡阶段；但"已修复"的记录是不准确的。

> **建议**：要么真正修 `common/scripts/run_md.py`（把 barostat 移到 `Simulation()` 之前，或加 `reinitializeContext()`），要么在 response/log 中把这条从"已修"改成"exp-A 专用脚本未修、影响已评估为可忽略"。**不要让一条未落地的修复留在"agreed/fixed"状态**——这正是 best-practice §22 警告的"SUPERSEDED 不标注会误导下一个读者"的变体。

**其它"待修"清单中仍未动的项（核实结果）：**
- `DPP4_FREE_RESIDUES` 集合仍有重复（`run_md.py:36`，631/547 各两次）——两份 review 都点过，仍在。
- `common/lib/` 仍为空，roadmap Phase 0.6 承诺的 `pdb_validator.py` 至今未实现——这是被自己写进路线图、又被 best-practice §6 强调的基础设施。
- 硬编码绝对路径、`paths.py` 未引入（`common/scripts/` 下只有 `build_system.py` 和 `run_md.py`）。

### 3.3 🔴 P0 安全问题以新形态复发：`cc.sh`

- May26 和 Jun01 两次点名 `cc-ds.sh` 明文 API key。核实：`cc-ds.sh` **已被 gitignore 且不在 git 历史中**（`git log --all -- cc-ds.sh` 为空）——这部分处理是到位的。
- **但 2026-06-06 新出现了 `cc.sh`**，含 token 字段，且 `git check-ignore cc.sh` 显示**它没有被 gitignore**，`git status` 里是 `?? cc.sh`（未跟踪但可被误 `git add`）。这是同一个 P0 问题换了个文件名复发。

> **建议**：(1) 立即把 `cc.sh` 加入 `.gitignore`（或并入已被忽略的模式）；(2) 确认其中 token 未被 `git add`；(3) 根治方案——把所有此类启动脚本统一改为从 `.env`/环境变量读取，杜绝"换个文件名又泄露"。

### 3.4 🟡 文档与实际状态的系统性漂移

随着项目推进，多处文档与实际状态脱节，建议做一次"文档-数据一致性审计"：

- **roadmap 进度表停留在 2026-05-27**：`roadmap.md:270` 标"最后更新 2026-05-27"，把 exp-C/exp-D 都标成"✓ 完成"，但其下"待启动"一行仍写"exp-D (Linker), exp-E, exp-F"——exp-D 同时是"完成"和"待启动"，自相矛盾。
- **exp-C README 的"成功判据"未随假说证伪而更新**：`exp-C/README.md:90` 仍写"C18 linker-二酸实现双点静电锚定"为成功判据，而实验结论恰恰是"linker 接上 FA 就逃逸、双点锚定不成立"。Jun01-response 已把"正式修订假说"列为 P1 agreed，但 README 正文未改。读者只看 README 会得到与实际相反的预期。
- **exp-A README 力场写 ff19SB，实际用 ff14SB**：`exp-A/README.md:57` MD 方案表写"ff19SB"，但同目录 `exp-log.md:18` 明确记录已改 ff14SB（理由：ff19SB CMAP 不支持 Aib）。决策是对的，README 表格没同步。
- **exp-A tasks.md 全部是 ⬜ 未做**：而 exp-log 显示 200 ns 生产 + 全量分析早已完成。tasks.md 完全没勾，已成僵尸清单。
- **exp-E / exp-F 目录为空**：roadmap 把它们列为项目"顶峰实验"，但无任何 README/占位。当前阶段可接受，但建议至少放一个 `README.md` 说明"未启动 + 前置依赖"。

### 3.5 🟡 跨实验方法学一致性需要在文档里讲清楚

- **系综不一致**：exp-A 用 NVT production，exp-C/D 用 NPT production。Jun01-response §3.4.3 说这是有意为之且"documented"，但我在各 README 里没找到集中、显式的说明。一旦未来做跨实验比较（尤其 exp-E 的 SAR 复现要把多个体系的能量放一起），系综差异必须在 Methods 里交代。
- **exp-D 体系无 HSA、exp-C 体系无受体**：这是合理的"分而治之"，但意味着**项目目前没有任何一个体系跑的是完整司美格鲁肽**（exp-C log 明确记"full sema restraint-based MD 16 轮全 NaN，不可行"）。这是一个需要在 roadmap 里诚实标注的能力边界——四步决策链里"脂链 vs 受体的空间竞争"这一核心矛盾，只有在完整分子（或 exp-F 的膜体系）里才能真正观察到。

### 3.6 🟡 分析统计严谨性（best-practice 自己的铁律仍未执行）

两份 Kimi review 都指出，至今仍未落实：
- `dt_ns` 硬编码而非从 reporter 间隔推导；
- 无自相关校正 / 有效样本量（best-practice §32 明文禁止对自相关序列直接做 mean±std 当作可比量）；
- 距离分布很可能多峰/偏斜，仍只报 mean±std，未给中位数/IQR 或核密度图；
- RMSD 未与 `gmx rms`/MDAnalysis 独立核对（best-practice §16）。

这些在 exp-D 的结论强度问题里直接放大了——**正因为没有统计校正，0.x Å 的差异才会被误读为"趋势"。** 这是优先级应当提到 P1 的项，而非一直 P2。

---

## 四、优先级建议

| 优先级 | 问题 | 建议 |
|--------|------|------|
| **P0（立即）** | `cc.sh` 含 token 未被 gitignore | 加入 .gitignore；确认未 `git add`；迁移到 `.env`/环境变量 |
| **P0（立即）** | exp-D 结论被高估（无统计、无 MM-GBSA、实际只跑 ~100ns） | 把 commit/README/log 的"qualitative match"加"not statistically significant"限定；规划补全 replica + correlated t-test + MM-GBSA |
| **P1（本周）** | barostat 时序"已修"但 `common/scripts/run_md.py` 实际未改 | 真正修复，或把记录从"已修"降级为"未修+影响可忽略"并标注 |
| **P1（本周）** | exp-C README 成功判据与证伪结论相反 | 正式修订假说，把 linker-FA 逃逸写成新发现（兑现 Jun01-response 的 agreed 项）|
| **P1（本周）** | 文档-数据一致性漂移（roadmap 进度、ff14/19SB、tasks 僵尸、500ns vs 100ns） | 做一次集中审计，统一 README/roadmap/tasks/exp-log |
| **P1（本周）** | 分析统计未达 best-practice §16/§32 | 自相关校正 + replica CV + 中位数/IQR + 独立工具核对 RMSD |
| **P2（生产前）** | `common/lib/pdb_validator.py` 未实现、paths.py 未引入、`DPP4_FREE_RESIDUES` 重复、旧脚本未归档 | 按既有清单逐项清 |
| **P2（持续）** | 跨实验系综/体系边界未集中说明 | 在 roadmap 或一份 Methods 草稿里集中交代 NVT/NPT、各体系缺失的组分、无完整 sema 体系的能力边界 |

---

## 五、结语

这个项目最珍贵的东西——清晰的科学叙事、诚实的失败记录、成熟的方法论自觉——依然完好，而且在 exp-B 暂停、exp-C 证伪、exp-D 的 NZ-C 键修复这几件事上得到了真实检验。这是少有的、值得长期信任的项目底色。

但当前阶段的主要矛盾，是**"推进速度"领先于"结论沉淀"**：exp-D 在几何指标差异不显著、结合自由能未算、模拟量只有计划 1/5 的情况下，就被写成了"复现 Novo Nordisk 设计逻辑"的成品 commit；同时，历史 review 里"已修"的 barostat 和反复点名的 API key，一个没真改、一个换名复发。这两类问题的共同根因是同一个：**记录的状态跑在了实际状态前面。** 而项目自己的 best-practice §22/§40 恰恰是为这件事写的——"诚实记录一切""被推翻的结论要标注，不要让下一个读者被误导"。

最务实的一步，不是再开 exp-E/F，而是回头把 exp-A–D 的结论强度、文档一致性、和那条没落地的修复**对齐到数据真实状态**。把 exp-D 的"qualitative match"诚实降级，把 exp-C 的假说正式修订，把 `cc.sh` 关掉——这三件事做完，项目的可信度会比再多跑几微秒 MD 提升得多。

---

*Review 完成时间：2026-06-06*
*Reviewer：Claude（Opus 4.8）*
*核对方式：18 个 md 文件全文 + 对 barostat/API key/力场参数/exp-D 数据的源码与文件系统交叉核对*
