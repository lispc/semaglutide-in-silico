# 司美格鲁肽 in silico 项目 Review — Kimi, 2026-05-26

> 本 review 基于对项目所有 Markdown 文档、Python/Shell 脚本、Amber/GROMACS 参数文件及实验记录的完整阅读。涵盖 `roadmap.md`、`best-practice.md`、实验 A 的全部文档与代码、通用脚本、参数文件等。

---

## 一、总体评价

这是一个**科学动机清晰、文献调研深入、文档质量极高**的计算结构生物学项目。作者明确知道自己的科学问题是什么（复现 Novo Nordisk 将 GLP-1 改造为司美格鲁肽的四步决策链），也清楚自己的贡献边界在哪里（已有文献各自覆盖片段，但缺乏统一的计算叙事）。从项目管理、实验设计到经验教训总结，都体现出前两个项目（cGAS-TRIM41、CJC-1295）积累的成熟方法论。

然而，项目在**代码安全性、脚本健壮性、力场参数完整性、GROMACS 协议细节、版本管理**等方面存在明显问题，部分问题可能在后续大规模生产阶段造成系统性错误甚至数据泄露。

---

## 二、做得好的地方 ✅

### 2.1 文档与知识管理（项目标杆级）

- **`roadmap.md` 逻辑链条极为清晰**：将司美格鲁肽的设计拆解为环环相扣的四步推理（Aib8 → K34R → C18 二酸 → γGlu-2×OEG），每一步都对应明确的计算实验（A–F）。这种"决策链驱动"的实验设计在学术计算项目中非常少见。
- **`best-practice.md` 是宝贵的组织资产**：基于约 5.5 μs 真实模拟经验总结的 40 条规范，涵盖引擎选择、力场处理、跨引擎验证、统计严谨性、错误检测、文档规范等。这不是泛泛而谈，而是带有具体数值教训的实战手册（如 `nstlist=40` 带来 23× 性能提升、RMSD 被高估 82% 的 Python 对齐 bug 等）。
- **`exp-log.md` 遵循"只追加、不删除"原则**：按时间顺序记录，每次记录含日期、操作、定量结果和决策依据。22 ns 的初步分析就给出了清晰的对比表格（WT vs Aib8 催化距离 5.2 vs 7.2 Å），并正确解读为"方向正确"。
- **中英文文档并存**：`docs/司美格鲁肽药物发现模拟设计.md` 是一份高质量的科普/综述性文档，包含详尽的 SAR 数据表和四个 in silico 实验方案，适合对外展示或论文 Introduction 的基础。

### 2.2 实验设计的科学严谨性

- **假说-验证结构明确**：每个实验都有核心假说、预期结果、成功判据和资源估算。例如实验 A 的成功判据包括 "Ser630 OG 到 Ala8 羰基碳平均距离 ≤ 3.5 Å"、"ΔΔG > 0" 等，可量化、可证伪。
- **Replica 设计不被妥协**：`roadmap` 和 `tasks.md` 中每个条件均规划 ≥3 replica，并明确要求报告 replica 间变异系数。
- **跨引擎验证（OpenMM + GROMACS）**：实验 A 同时启动了 OpenMM（主引擎，有约束）和 GROMACS（验证引擎，无约束）的模拟。这直接体现了 `best-practice.md` §7 的精神。
- **渐进式构建**：先跑 200 ns × 1 replica 验证体系稳定性，确认方向正确后再投入 3×300–500 ns 的生产。这避免了"一步到位后发现体系崩溃"的灾难。
- **约束策略经过思考**：针对 CJC-1295 项目"长肽从口袋漂移"的教训（Trap 5/6），设计了分区域约束——DPP-4 CA 弱约束 + 肽 N 端 backbone 中等约束 + C 端自由 + 催化位点完全不约束。这在保持关键区域自由度的同时解决了漂移问题。

### 2.3 代码的关键亮点

- **`run_md.py` 的协议完整**：包含 0→100 K NVT 加热、100→310 K NPT 渐进加热、NPT equilibration、NVT production 的完整四阶段协议；支持 `--restart` 从 checkpoint 续跑；有进度报告（ns/day、剩余时间）。
- **`build_system.py` 的错误恢复机制**：对 `modeller.addHydrogens()` 的 `ValueError` 做了 retry 逻辑（删除模板不匹配的残基后重试），这是处理非标准 termini 的实用技巧。
- **`quick_check.py` 的鲁棒性设计**：当 `resid 630 and name OG` 选择失败时（tleap 可能改变编号），自动 fallback 到"找离肽最近的 SER OG"的策略，体现了对 PDB 编号不可靠性的防御式编程。
- **`setup_gmx.py` 的 MDP 自动生成**：将 Amber 体系转换为 GROMACS 时自动写入 em/nvt/npt/md 四组 mdp 文件，参数与 OpenMM 协议对应，降低了手动编辑出错概率。

### 2.4 项目管理与资源规划

- **清晰的目录结构**：`exps/exp-A/` 下 `structures/`、`tleap/`、`md/`、`gmx/`、`analysis/` 分工明确；`common/params/` 和 `common/scripts/` 实现了跨实验共享。
- **资源估算现实**：`roadmap.md` 中给出了每个实验的体系数、单体系模拟量、replica 数、总模拟量和日历天数（按 4×RTX 3090 估算），总计约 47 μs / 35 天，对单课题组项目而言是合理且可执行的。
- **`.gitignore` 设计合理**：排除了轨迹（DCD/XTC）、checkpoint、大型拓扑文件、GROMACS 输出等，同时保留了代码和文档。

### 2.5 对文献的精准定位

- `roadmap.md` 中对已有文献（Frimann 2023、Liu 2025、Sønderby 2016 等）的"做了什么 / 没做什么"表格，精准地划出了本项目的学术空间。这种诚实的文献定位避免了"重复造轮子"的同时，也明确了真正的创新点。

---

## 三、做得不好甚至不对的地方 ❌

### 3.1 🔴 安全与隐私（严重）

- **`cc-ds.sh` 硬编码 API Key**：
  ```bash
  export ANTHROPIC_AUTH_TOKEN=sk-27ce3e9b64794ee886b7c1077a44e8a2
  ```
  这是一个明文写入仓库的 DeepSeek API Key。尽管 `.gitignore` 排除了 `cc-ds.sh`，但一旦曾经 `git add` 过，该 key 就永久留在 git 历史中（`git log -p` 仍可提取）。**此 key 应立即在 DeepSeek 控制台撤销并更换**，且未来所有 API key 应通过环境变量或 `.env` 文件（加入 `.gitignore`）注入。

### 3.2 🔴 力场参数不完整（可能影响模拟正确性）

- **`common/params/aib_residue.xml` 只有 Bond，没有 Angle/Dihedral/Improper**：
  ```xml
  <Bond atomName1="CA" atomName2="CB"/>
  <Bond atomName1="CA" atomName2="CB2"/>
  <!-- 无 <Angle>、<Proper>、<Improper> 定义 -->
  ```
  OpenMM 的 `ForceField` 对非标准残基如果模板中只定义了 Bond，会从力场中匹配 generic angle/dihedral 参数。但 Aib（α-氨基异丁酸）的侧链是两个甲基连在同一个 CA 上，其 CA-CB 键长、CA-CB-HB 角度、以及 backbone 附近的二面角分布与标准氨基酸有显著差异。**缺少明确的 angle 和 dihedral 定义可能导致 Aib8 区域的局部几何弛豫到非物理状态。**
  - 建议：补充从 GAFF2 (`aib.frcmod`) 中提取的 angle、proper dihedral 和 improper dihedral 参数到 XML 中；或在 `build_system.py` 中验证 Aib8 区域的键长/角度分布是否合理。

- **Aib 的电荷来源说明与实现不完全一致**：
  `aib_residue.xml` 的注释说 "charges from alanine"，但 Ala 的总电荷为 0（N: -0.4157, H: 0.2719, CA: 0.0337, HA: 0.0823, CB: -0.1825, HB1/2/3: 0.0603, C: 0.5973, O: -0.5679）。Aib 去掉了 HA（+0.0823），增加了 CB2（-0.1721）和 HB21/22/23（3×0.0707）。XML 中的 CA 电荷改为 0.0341（与注释一致），但 CB 的电荷从 -0.1825 改为 -0.1721，这实际上是**从 Ala 的电荷重新分配而来**，而非直接使用。总和约为 0，但需要明确说明电荷重分配的逻辑。

### 3.3 🟡 GROMACS 协议缺陷

- **缺少渐进加热**：
  `setup_gmx.py` 生成的 `nvt.mdp` 在 100 K 跑 50 ps，`npt.mdp` 直接跳到 310 K 跑 300 ps。中间没有 100→310 K 的渐进升温步骤。虽然 300 ps NPT 后体系可能已经平衡，但**直接从 100 K 跳到 310 K 对蛋白体系是危险的**，可能导致局部热冲击或初期结构畸变。
  - 建议：NVT 阶段做分步升温（100→150→200→250→310 K），或至少用 `gen_vel` 重新生成速度后短平衡。

- **`emtol = 1000` 过松**：
  `mdp_em` 中 `emtol = 1000`（kJ/mol/nm）在能量最小化中过于宽松。标准做法通常是 `emtol ≤ 100`，甚至 ≤ 10，以确保初始构象没有严重的局部应力。1000 的阈值可能掩盖 frame-0 的 clash。

- **NPT equilibration 仅 300 ps**：
  对于 ~110–140k atoms 的蛋白-水体系，300 ps 的 NPT 平衡通常不足以让密度和盒体积完全收敛。建议至少 1–2 ns。

- **缺少 `lincs_iter` / `lincs_order` 优化**：
  `best-practice.md` §8 明确建议 `lincs_iter=2, lincs_order=6`，但生成的 mdp 文件未设置这些参数（使用 GROMACS 默认值 `lincs_iter=1, lincs_order=4`）。对于 200+ ns 的长模拟，默认约束精度可能不足。

- **GROMACS 交叉验证"无约束"导致不公平比较**：
  `exp-log.md` 明确记录 "GROMACS: 无约束（仅用作引擎交叉验证）"。但 OpenMM 施加了 DPP-4 CA 弱约束 + 肽 N 端 backbone 约束。**两个引擎的体系约束条件不同，使得"跨引擎验证"的比较基础不公平**。验证的目的应是确认"同一物理体系在不同数值实现下行为一致"，而约束差异引入了额外的物理差异。
  - 建议：GROMACS 交叉验证也应施加对应的 position restraints（`posres.itp`），或在分析时明确说明约束差异对比较的影响。

### 3.4 🟡 代码健壮性与工程债务

- **多处硬编码绝对路径**：
  `common/scripts/build_system.py`、`run_md.py`、`exps/exp-A/analysis/quick_check.py`、`exps/exp-A/tleap/` 下的所有脚本，均硬编码了 `/home/scroll/personal/semaglutide-in-silico`。这使得代码无法在其他机器或用户目录下运行，也阻碍了 `common/scripts/` 被 exp-B–F 复用。
  - 建议：使用 `pathlib.Path(__file__).resolve().parents[N]` 或环境变量推导项目根目录。

- **`run_md.py` 硬编码 `CudaDeviceIndex: '0'`**：
  脚本不支持指定 GPU ID，多任务并行时必须手动修改源码或复制脚本。`run_gmx.sh` 通过 `$GPU` 参数支持多 GPU，但 OpenMM 端缺失对应能力。

- **`build_system.py` 第 63 行有未完成的 `pass` 分支**：
  ```python
  if not to_delete:
      for chain in modeller.topology.chains():
          residues = list(chain.residues())
          pass  # <-- 空逻辑
      raise
  ```
  这段代码是在处理 `addHydrogens` 失败时尝试删除最后一个残基，但 `pass` 分支意味着如果按 `res.index` 找不到残基，就直接 `raise`，而不尝试备用策略。这是一个未完成的错误恢复逻辑。

- **`DPP4_FREE_RESIDUES` 集合中有重复值**：
  ```python
  DPP4_FREE_RESIDUES = {630, 708, 740, 547, 629, 631, 656, 662, 666, 710, 711, 631, 632, 633, 634, 635, 547, 548}
  ```
  `631` 和 `547` 各出现了两次。集合会自动去重，但这反映出代码审查不够仔细。

- **三个不同版本的 Aib8 修改脚本共存，未归档**：
  `tleap/` 下同时存在 `convert_to_aib8.py`、`fix_aib8_prmtop.py`、`fix_aib8_v2.py`。根据 `exp-log.md` 和 `minimize_aib8.py` 的引用，实际使用的是 `fix_aib8_prmtop.py`（生成 `aib8_modified.prmtop`）。旧版本应按 `best-practice.md` §18 的建议归档到 `archive/_versions/`，而不是留在生产目录中造成混淆。

- **`setup_gmx.py` 使用 parmed 转换 Amber→GROMACS**：
  `best-practice.md` §8 明确警告 "parmed 的 AMBER→GROMACS 转换丢失了残基特异性 CMAP"，并建议 "始终用原生实现作为 gold standard"。虽然目前使用 ff14SB（无 CMAP），所以当前无影响，但如果后续实验（如实验 F 的全长 GLP-1R-TMD）改用 ff19SB，这个转换策略将引入系统性误差。**这是一个被自己的 best-practice 标记为陷阱的做法。**

### 3.5 🟡 实验设计与文档一致性问题

- **模拟时长在多个文档中不一致**：
  - `roadmap.md`：实验 A 为 "500 ns cMD × 3 replica × 2 systems = 3 μs"
  - `exp-A/README.md`：第一阶段 200 ns × 1 replica，第二阶段 300–500 ns × 3 replica
  - `exp-A/tasks.md`：Phase 2 是 "各 200 ns, 1 replica"，Phase 5 是 "× 2 additional replica × 300 ns"
  实际总模拟量（如果执行 Phase 5）约为 2×(200 + 3×300) = 2.2 μs 或 2×(200 + 3×500) = 3.4 μs，与 roadmap 的 "3 μs" 大致吻合，但各文档的表述不一致可能导致执行困惑。

- **`tasks.md` 中任务编号重复**：
  存在两个 "A-23"（A-23 肽 backbone RMSF 对比，A-23b C 端漂移检查）。编号系统应唯一。

- **`roadmap.md` Phase 0.6 承诺的 `lib/pdb_validator.py` 未实现**：
  `common/lib/` 目录为空，但 `best-practice.md` §6 和 `roadmap.md` 均强调需要一个经过列验证的 PDB writer/validator。这是一个被记录为教训但尚未修复的问题。

- **OpenMM 的 NPT equilibration 仅 200 ps**：
  `run_md.py` 中 NPT 阶段为 `simulation.step(100000)` = 200 ps。与 GROMACS 端类似，这个时长对 140k atoms 体系偏短。建议在 200 ps 后检查密度/体积收敛曲线。

### 3.6 🟡 分析脚本的统计与方法论问题

- **`quick_check.py` 硬编码 `dt_ns = 0.1`**：
  这个值应该从 DCD reporter 的间隔（`dcd_interval = 50000` steps × 2 fs = 100 ps = 0.1 ns）推导，而不是硬编码。如果将来改变输出频率，分析脚本会默默产生错误的时间轴。

- **未做自相关校正**：
  `quick_check.py` 报告 catalytic distance 的 "mean ± std"，但 MD 轨迹是高度自相关的时间序列。`best-practice.md` §32 明确警告 "不能用独立样本 t-test"，但分析脚本目前只做了最简单的均值/标准差计算，没有有效样本量校正或 blocked bootstrap。

- **距离分布的描述统计选择**：
  催化距离（Ser630 OG → Ala8/Aib8 C）的分布在 22 ns 数据中范围为 3.7–8.2 Å（WT）和 5.4–8.6 Å（Aib8），这种分布很可能是多峰或严重偏斜的。报告 "mean ± std" 可能误导读者假设其服从高斯分布。建议同时报告中位数、IQR，或画完整的核密度估计图。

- **RMSD 计算未验证对齐方法**：
  `quick_check.py` 使用 `mdtraj` 的 `superpose()` 对齐，这与 `best-practice.md` §16 的 "分析脚本必须用独立方法验证" 精神相悖。虽然 `mdtraj` 是成熟库，但在正式结论前应与 `gmx rms` 或 MDAnalysis 的 `alignto()` 核对一次。

### 3.7 🟡 结构文件与编号管理

- **大量中间 PDB 文件缺乏版本控制**：
  `structures/` 下有 `DPP4_GLP1_Aib8.pdb`、`_v2.pdb`、`_noh.pdb`、`_clean.pdb`、`_final.pdb`、`_tleap.pdb`、`_hybrid.pdb`、`_withH.pdb`、`_asALA_noh.pdb`、`DiprotinA_ref.pdb` 等 20+ 个文件，但没有文档说明每个文件的生成顺序、用途和 "gold standard" 是哪个。这与 `best-practice.md` §15 "建立编号字典" 的建议相矛盾。

- **`quick_check.py` 硬编码 resid 编号**：
  ```python
  dpp4_idx = t.topology.select("resid 0 to 727 and backbone")
  pep_idx = t.topology.select("resid 728 to 758 and backbone")
  ```
  这些编号依赖于 tleap 加载 chain B + peptide 的具体顺序。如果未来构建其他体系（如实验 B 的 GLP-1R ECD），这些硬编码编号将失效。建议用 `chain.id` 或残基名选择，而非绝对 resid。

### 3.8 🟡 力场选择记录的内部不一致

- `roadmap.md` 技术规范表中写 "蛋白力场: Amber14SB **或 ff19SB**"，但 `exp-A/README.md` 和 `exp-log.md` 明确使用 ff14SB，理由是 "ff19SB 的残基特异性 CMAP 类型不支持非标准 Aib 残基"。这个决策是合理的，但 `roadmap.md` 的 "或 ff19SB" 表述可能在后续实验中造成混淆——如果实验 F 使用全长 GLP-1R-TMD 并考虑 ff19SB，需要重新评估 Aib 的兼容性。

---

## 四、优先级修复建议

| 优先级 | 问题 | 建议修复 |
|--------|------|----------|
| **P0 (立即)** | API Key 泄露 (`cc-ds.sh`) | 1. 在 DeepSeek 控制台撤销该 key；2. 使用 `.env` 文件 + `python-dotenv` 或 shell 环境变量注入；3. 考虑 `git filter-branch` 或 BFG Repo-Cleaner 从历史中移除该文件 |
| **P0 (立即)** | `aib_residue.xml` 缺少 Angle/Dihedral | 从 `aib.frcmod` / `aib_capped.frcmod` 提取 angle、proper dihedral、improper dihedral 参数，补充到 XML 模板中 |
| **P1 (本周)** | GROMACS 加热协议缺陷 | 在 `setup_gmx.py` 的 `nvt.mdp` 中实现分步升温，或增加 `annealing` 协议；`emtol` 降至 ≤100 |
| **P1 (本周)** | 硬编码绝对路径 | 统一改为 `pathlib` 相对路径或环境变量；`common/scripts/` 应通过 `argparse` 接收 `--exp-dir` 参数以支持 exp-B–F 复用 |
| **P1 (本周)** | GROMACS 交叉验证约束不一致 | 为 GROMACS 生成对应的 `posres.itp` 文件，使约束条件与 OpenMM 一致 |
| **P2 (生产前)** | `quick_check.py` 统计方法升级 | 增加有效样本量校正、非高斯分布的中位数/IQR 报告、自相关时间估计 |
| **P2 (生产前)** | 清理/归档旧版脚本 | 将 `convert_to_aib8.py`、`fix_aib8_v2.py` 移入 `archive/`；在 `README.md` 中说明当前使用的构建流程 |
| **P2 (生产前)** | 实现 `pdb_validator.py` | 按照 `best-practice.md` §6 的建议，写一个列对齐验证器，统一所有脚本的 PDB I/O |
| **P3 (持续)** | 文档一致性审计 | 统一 `roadmap.md`、`README.md`、`tasks.md` 中的模拟时长、replica 数、力场选择等关键参数 |

---

## 五、结语

这是一个**起点很高、潜力巨大**的项目。从科学问题到实验设计，从文献调研到经验总结，都展现出成熟研究者的方法论水平。当前的主要风险集中在**代码安全（API key）、力场参数完整性和协议细节**三个层面——前两者是"做对了就锦上添花，做错了就全盘皆输"的类型。

特别值得肯定的是 `best-practice.md` 的建立。如果项目团队能够**严格遵循自己写下的规范**（如 §16 的分析脚本验证、§8 的 parmed 警告、§35 的 frame-0 检查清单），并及时修复当前发现的偏差，那么后续 Phase 1–2 的大规模生产将建立在非常坚实的基础上。

---

*Review 完成时间：2026-05-26*
*Reviewer：Kimi*
