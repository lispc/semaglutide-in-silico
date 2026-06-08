# 司美格鲁肽 in silico 项目 Review — Kimi, 2026-06-01

> 本 review 基于对项目所有 Markdown 文档、Python/Shell 脚本、Amber/GROMACS 参数文件及实验记录的完整阅读。涵盖 `roadmap.md`、`best-practice.md`、实验 A–D 的全部文档与代码、通用脚本、参数文件等。与 2026-05-26 的 review 相比，项目已完成了 exp-A 的全部 200 ns 生产模拟，exp-C 的 9 条 100 ns 游离 FA/linker 验证，以及 exp-D 的 5 种 linker 变体构建与 MD 运行。

---

## 一、总体评价

项目在**过去 5 天内取得了显著的实验进展**：exp-A 的 200 ns×2 体系生产完成并给出了方向正确的核心结论；exp-C 的游离 FA 实验揭示了 unexpected 的 linker-FA 逃逸现象，这是一个具有科学价值的负面结果；exp-D 的 5 种 linker 变体已启动 MD 生产。然而，**代码质量、工程债务和文档一致性在快速推进中出现了明显恶化**，特别是在 exp-D 中出现了多个可能被判定为"会直接导致错误数据"的 bug。

当前项目处于"高速推进但地基在摇晃"的状态。如果不对 exp-D 的构建脚本和 launcher 进行紧急修复，后续分析可能基于物理上错误的拓扑文件或无法运行的模拟。

---

## 二、做得好的地方 ✅

### 2.1 科学决策的理性与灵活（项目级亮点）

- **exp-B 的暂停决策堪称典范**：在发现所有已发表晶体结构（3IOL/4ZGM/6GB1/7KI0）中 ECD 与肽链均分离 25–40 Å 后，团队没有强行推进 FlexPepDock（已在 10 decoy 中验证无效），而是果断暂停 exp-B、转向 exp-C。这体现了 "知道什么时候停" 的 best-practice §29 精神。
- **exp-C 的迭代式参数化展示了优秀的 debug 能力**： linker-C18 经历了四轮电荷方案迭代（全局 charge scaling → GAFF2 逐原子 → RDKit extended conformer + 分级加热 → AM1-BCC），最终定位到"手工 GAFF2 电荷过于极化"的根因。这种系统性排除法值得肯定。
- **exp-C 的负面结果被如实记录**："只要 linker 接上，FA 就呆不住 FA3"这一与原始假说相反的结论，在 exp-log.md 中被完整记录而非掩盖。这符合 best-practice §40 "诚实记录一切" 的原则。

### 2.2 实验进展与数据产出

- **exp-A 完成 200 ns 生产**：OpenMM WT/Aib8 各 200 ns + GROMACS 交叉验证各 200 ns。核心结论清晰：Aib8 双甲基在所有 5 个指标上推离 DPP-4（催化距离 +1.0 Å, 接触数 -49, S1 距离 +1.7 Å）。
- **exp-C 完成 9 条 100 ns 轨迹**：游离 C18 单酸×3、游离 C18 二酸×3、linker-C18×3。发现远端羧基在游离 FA 中稳定锚定 ARG482（2.8 Å），但在 linker 接上后逃逸至 32–41 Å。
- **exp-D 启动 5 种 linker 变体 MD**：no_linker、γGlu、γGlu-1×OEG、γGlu-2×OEG、γGlu-3×OEG，每种 3 replica，预计总模拟量 ~7.5 μs。

### 2.3 文档与知识管理（总体维持高标准）

- **`roadmap.md` 及时更新进度**：2026-05-27 的"项目当前进度"表格清晰标注了 exp-A(完成)、exp-C(完成)、exp-D(运行中)、exp-B(暂停) 的状态。
- **`exp-log.md` 依然是最佳实践范例**：exp-C 的日志按时间顺序记录了 NaN 诊断（4 轮迭代）、AM1-BCC 修复、ACE-linker 假说证伪等完整决策链。
- **`best-practice.md` 的 40 条规范仍在指导实践**：exp-A 中 `nstlist=40` 的性能修复、exp-B 中 FlexPepDock 不可靠的判断、exp-C 中 frame-0 clash check 等，均能看到 best-practice 的影子。

### 2.4 代码中的关键亮点

- **`run_md.py`（exp-C/exp-D 版）的通用性提升**：相比 exp-A 的硬编码 WT/Aib8，exp-C/D 的 `run_md.py` 通过 `--system` 参数支持多体系复用，是一个明显的工程进步。
- **`build_lya_variants.py` 的批量参数化设计**：用 SMILES 片段 + RDKit ETKDGv3 + antechamber AM1-BCC 的流水线，一次性生成 5 种 linker 变体的 mol2，设计思路清晰。
- **`analyze_ace_linker.py` 的鲁棒性**：使用 `md.iterload` 处理不完整的 DCD，对 rep1 DCD 损坏的情况 gracefully 跳过。

---

## 三、做得不好甚至不对的地方 ❌

### 3.1 🔴 会直接导致错误数据/模拟失败的 bug

#### 3.1.1 `exps/exp-D/md/launch.sh`：跨实验路径错误 + CLI 参数不匹配

```bash
# 第 38 行
python -u "$REPO/exps/exp-C/md/run_md.py" \
    --system "$vname" --replica $rep --gpu 0 \
    --prmtop "$TLEAP/$vname.prmtop" \
    --inpcrd "$TLEAP/$vname.inpcrd" \
    --outdir "$outdir"
```

- **路径错误**：exp-D 的 launcher 调用的是 **exp-C 的 `run_md.py`**，而非 exp-D 自己的 `run_md.py`。虽然两个文件内容相似，但 exp-C 的 `run_md.py` 硬编码了 `EXP_C` 路径，会导致所有输出被写入 `exps/exp-C/md/...` 而非 `exp-D/md/...`。
- **参数不匹配**：exp-C 的 `run_md.py` 只接受 `--system`、`--replica`、`--gpu`、`--nsteps` 四个参数。`--prmtop`、`--inpcrd`、`--outdir` 是不存在的参数，会导致 `argparse` 报错并立即退出。
- **GPU 分配错误**：`CUDA_VISIBLE_DEVICES=$gpu` 设置了正确的 GPU，但 `--gpu 0` 又把 OpenMM 的 CudaDeviceIndex 固定为 0。如果 `CUDA_VISIBLE_DEVICES=1`，OpenMM 会看到 "唯一的 GPU 是 index 0"，这恰好能工作，但如果 `CUDA_VISIBLE_DEVICES` 被省略或脚本在其他机器上运行，会导致所有任务挤在 GPU 0 上。

**后果**：此 launcher 无法成功启动任何 exp-D 的 MD 模拟。脚本会立即因未知参数而退出。

#### 3.1.2 `exps/exp-D/build/build_complex.py`：Lys→Arg 突变逻辑化学不正确

```python
if a.name == 'CE':
    a.name = 'NE'
    a.type = 'NE'
elif a.name == 'NZ':
    a.name = 'CZ'
    a.type = 'CZ'
```

- **问题**：LYS 的 CE（ε-碳）和 NZ（ζ-氮）与 ARG 的 NE（ε-氮）和 CZ（ζ-碳）在化学上是**不同的元素和几何位置**。简单重命名原子标签而不改变坐标，意味着 ARG 的胍基被放在了 Lys 侧链的线性胺位置上——几何完全错误（ARG 胍基是平面三取代，LYS 胺是四面体）。
- **tleap 后果**：tleap 加载这个 PDB 后，会根据 ARG 的模板重新构建缺失原子（NH1, NH2, NE, CZ 的 H）。但由于骨架原子（NE/CZ）的初始位置是错的，tleap 构建的新原子会基于错误的几何生成，导致局部结构畸变。

**后果**：exp-D 中所有体系的 Lys34→Arg 突变在拓扑层面是**化学上不正确的**。如果 Arg34–Glu27 的氢键网络是实验核心指标，这个错误的突变会使结论不可信。

#### 3.1.3 `exps/exp-D/build/build_complex.py` 生成的 tleap 脚本：双键 NZ

```tcl
bond complex.116.NZ LYA.1.N
bond complex.116.NZ LYA.1.C
```

- **问题**：同一个 NZ 原子被同时 bond 到 LYA 的 N 和 C 两个原子上。这在化学上意味着 NZ 形成了**三键（一个肽键 + 一个额外的 C-N 键）**，或者更准确地说是**四配位氮**——这在标准氨基酸化学中是不存在的。
- **tleap 后果**：tleap 可能会静默接受这个指令并生成一个奇怪的拓扑，或者因 valence 错误而失败。如果静默成功，模拟中 NZ 的几何会因错误的键连接而严重畸变。

**后果**：exp-D 的 `gglu_1oeg`、`gglu_2oeg`、`gglu_3oeg` 变体（使用 `build_gglu_1/2/3oeg.in`）的拓扑文件可能包含**化学上无意义的多配位氮**。这可能导致能量爆炸、NaN 或虚假的构象约束。

#### 3.1.4 `common/scripts/run_md.py`：Barostat 在 Simulation 创建后才添加

```python
simulation = app.Simulation(...)
# ... 若干行后 ...
system.addForce(mm.MonteCarloBarostat(...))
```

- **问题**：OpenMM 的 `System` 在传入 `Simulation` 构造函数后，**后续对 `System.addForce()` 的调用不会自动应用到已有的 `Context`**。`Simulation` 不会自动调用 `reinitializeContext()`。
- **后果**：NPT heating 和 NPT equilibration 阶段实际上运行在 **NVT（无 barostat）** 下。虽然 200 ps 的 NPT equilibration 可能因其他原因（如初始密度已合理）而看起来"正常"，但压力耦合实际上没有生效。这对于需要严格 NPT 收敛的体系是一个系统性缺陷。

#### 3.1.5 `exps/exp-A/gmx/run_gmx.sh`：坐标文件名不匹配

```bash
COORD="${LABEL}_coord.gro"
```

但 `setup_gmx.py` 生成的是 `${LABEL}.gro`。脚本会在 `gmx grompp` 阶段因找不到文件而失败，除非用户手动重命名。

### 3.2 🔴 力场参数与化学完整性问题

#### 3.2.1 `common/params/aib_capped.mol2`：化学错误

- **CH3–O 键形成三元环**：ACE cap 中多了一个不应该存在的 CH3-O 键，创建了一个化学上不可能的三元环结构。
- **C=O 键被标记为单键**：羰基双键类型错误为单键。
- **后果**：`aib_capped.frcmod` 从这个错误的 mol2 导出，包含物理上不可能的键角（如 `c2-c1-os = 179.5°`，在甲基碳上呈线性）。这些文件虽然当前未被直接使用，但存在于 `common/params/` 中，存在被误用的风险。

#### 3.2.2 `common/params/sqm.out`：AM1 优化失败

- 最终热函 **+257 kcal/mol**（异常高）。
- CA–CB2 键长膨胀至 **2.61 Å**（标准 C-C 单键 ~1.54 Å）。
- NME 氢原子脱落至空间完全无关的位置。
- **后果**：`sqm.pdb` 包含断裂的几何结构，绝不能用于下游参数化。虽然项目实际使用的是初始几何而非 sqm 输出，但这个失败记录提示 antechamber/sqm 管线在该分子上不可靠。

#### 3.2.3 `common/params/aib_residue.xml`：仍缺少 Angle/Dihedral（May26 已指出，未修复）

- 只有 `<Bond>` 定义，无 `<Angle>`、`<Proper>`、`<Improper>`。
- OpenMM 会尝试从力场中匹配 generic angle/dihedral，但 Aib 的双甲基几何与标准氨基酸显著不同。
- **当前声称"无实际影响"的理由不成立**：虽然 exp-A 的 200 ns 模拟未出现明显 crash，但这不等于 Aib8 区域的局部几何是物理正确的。键长/角度的细微偏差可能在 MM-PBSA 能量分解中累积为系统误差。

#### 3.2.4 `build_lya_v2.py`：RDKit 原子顺序硬编码

```python
bb_n=3; bb_ca=4; bb_c=57; bb_o=58; bb_cb=5; nme_atoms={59,60}
```

- 这些索引完全依赖于特定 RDKit 版本对特定 SMILES 的 canonical atom ordering。
- 如果 RDKit 版本升级或 SMILES 解析器行为改变，所有 backbone 原子名称都会错配。
- **后果**：LYA 残基的 backbone 可能被错误地标记，导致 tleap 无法识别 N/Cα/C 位置，或产生错误的二面角约束。

#### 3.2.5 `build_sema_parmed.py`：通用键类型替代所有未定义键

```python
default_type = pmd.BondType(300.0, 1.40)
```

- 用一个通用的单键参数（k=300, r=1.40）替代 linker 中所有未被 antechamber/tleap 识别的键。
- linker 中包含酰胺键（C=O 双键特征，r~1.23 Å, k~600+）、醚键（C-O-C, r~1.43 Å）等多种键型。
- **后果**：模拟中 linker 的局部柔性可能被系统性地错误描述，影响脂链-肽构象动力学。

### 3.3 🟡 代码工程债务急剧累积

#### 3.3.1 硬编码绝对路径仍未解决

- 从 `common/scripts/` 到 `exps/exp-D/build/`，**所有 40+ 个脚本**仍硬编码 `/home/scroll/personal/semaglutide-in-silico`。
- `build_lya.py`、`build_lya_v2.py`、`build_combined_pdb.py` 等还硬编码 `os.chdir(...)`，修改全局工作目录。
- **后果**：代码完全不可移植；任何在新机器或不同用户目录下的复现尝试都会失败。

#### 3.3.2 exp-D 完全没有项目文档

- exp-D 目录下**不存在** `README.md`、`exp-log.md`、`tasks.md`。
- 这是与 exp-A/B/C 相比的**显著倒退**。exp-D 的科学目标（验证 γGlu-2×OEG 如何化解空间冲突）、实验设计、成功判据、资源估算均无从查证。
- 只能从 `build_lya_variants.py` 的注释和 `exp-C/exp-log.md` 的末尾推断 exp-D 的动机。
- **后果**：3 天后（甚至更短时间内），连作者自己都可能忘记 exp-D 中每个变体的具体构建流程和假设。

#### 3.3.3 `exps/exp-A/analysis/mmpbsa.py`：不完整脚本仍在生产目录

```python
# 文件末尾
print("\n=== Summary ===")
print(f"{'Metric':<30} {'WT':>12} {'Aib8':>12}")
# These will be filled after both loops
# PYEOF
```

- 汇总表格的标题被打印，但数据行缺失。脚本在运行到此行后结束，不输出任何定量比较结果。
- **后果**：如果有人运行此脚本期望获得 MM-PBSA 对比，会得到一个空表，误以为数据有问题。

#### 3.3.4 旧版本脚本未归档

- `exps/exp-A/tleap/` 下仍共存 `convert_to_aib8.py`、`fix_aib8_prmtop.py`、`fix_aib8_v2.py`、`minimize_aib8.py`。
- `exps/exp-C/tleap/` 下共存 `build_diacid.py`、`build_diacid_v2.py`、`build_diacid_pdb.py`、`build_combined_pdb.py` 等多个功能重叠的脚本。
- **后果**：生产目录混乱，无法确定哪个脚本是"gold standard"。

#### 3.3.5 大量 Python 脚本缺少 `if __name__ == '__main__'`

- `build_diacid.py`、`build_diacid_v2.py`、`build_combined_pdb.py`、`build_sema_pdb.py`、`build_sema_step2.py`、`fix_mol2.py`、`position_lnk.py`、`minimize_aib8.py` 等均在 import 时立即执行。
- **后果**：这些脚本不可被其他脚本安全地 import，也无法编写单元测试。

#### 3.3.6 `run_md.py`（common）中 `DPP4_FREE_RESIDUES` 集合重复

```python
DPP4_FREE_RESIDUES = {630, 708, 740, 547, 629, 631, 656, 662, 666, 710, 711, 631, 632, 633, 634, 635, 547, 548}
```

- `631` 和 `547` 各出现两次。集合会自动去重，但这反映出代码审查不够仔细。

### 3.4 🟡 实验设计与科学方法论问题

#### 3.4.1 exp-D 的体系构建策略存在根本性化学问题

- `build_complex.py` 试图通过 PDB 原子重命名实现 LYS→ARG 突变，但这不是一个有效的突变策略。
- 正确的做法应该是：使用 tleap 的 `mutate` 命令，或用 Modeller/PyMOL 做侧链替换后重新能量最小化。
- 当前方法可能导致 Arg34 的胍基处于错误的位置，从而使"Arg34–Glu27 氢键网络"这一核心分析指标失去物理意义。

#### 3.4.2 exp-C 的" linker-FA 逃逸"发现对原假说的冲击

- 原假说（roadmap.md、exp-C/README.md）：linker 是"主动定位器"，将近端羧基精确放置在 HSA 表面第二个正电荷残基附近，实现双点锚定。
- 实际发现：无论是否 ACE 封端，linker-FA 均从 FA3 逃逸（32–41 Å）。驱动力是 OEG 单元的亲水性，而非末端电荷。
- **问题**：这个负面结果本身是有价值的科学发现，但 roadmap.md 和 exp-C/README.md 中的"成功判据"（"C18 linker-二酸实现双点静电锚定"）可能需要被修订。目前没有看到对假说的正式修正。

#### 3.4.3 exp-D 的 MD 协议中未移除 Barostat

- `run_md.py`（exp-D 版）在 heating 后添加了 `MonteCarloBarostat`，然后在 production 阶段**没有移除**它。
- 虽然对于 NPT production 这是正确的（与 exp-A 的 NVT production 不同），但注释写的是 "NPT eq done, starting production"，未明确说明 production 是 NPT 还是 NVT。
- 更大的问题是：exp-A 声称使用 NVT production（见 `exp-A/README.md`），而 exp-D 使用 NPT production。不同实验使用不同系综，跨实验比较时需要谨慎。

#### 3.4.4 `addIonsRand Na+ 0 Cl- 0` 不加离子

- `build_complex.py`、`parmed_assemble.py`、`write_tleap.py` 及多个 `.in` 文件中均使用 `addIonsRand complex Na+ 0 Cl- 0`。
- 如果体系有净电荷（如 linker 的 -1 或 -2 羧酸根），此命令**不会添加任何抗衡离子**，导致模拟在净电荷非零的条件下运行。
- **后果**：未中和的体系在 PME 下会产生非物理的均匀背景电荷校正，可能影响静电相互作用。

#### 3.4.5 分析脚本的统计方法仍未升级

- `quick_check.py`、`full_analysis.py` 仍使用硬编码 `dt_ns = 0.1`，未从 DCD reporter 间隔推导。
- 未做自相关校正：报告 "mean ± std" 但未校正时间序列自相关。
- 未报告中位数/IQR：催化距离分布可能是多峰或偏斜的，仅用 mean±std 可能误导。
- RMSD 计算未与 `gmx rms` 或 MDAnalysis 独立验证。
- 这些在 May26 review 中已指出，至今未修复。

### 3.5 🟡 参数文件与配置问题

#### 3.5.1 `.claude/settings.local.json` 包含过宽的通配符和 stale PID

```json
{"command": "Bash(bash *)", "allowed": true},
{"command": "Bash(git *)", "allowed": true},
{"command": "Bash(python3 *)", "allowed": true},
{"command": "Bash(kill 3806129)", "allowed": true}
```

- `bash *` 和 `python3 *` 的通配符允许任意 shell/python 命令执行，实质上绕过了权限白名单。
- `kill 3806129` 是一个已不存在的硬编码 PID，是 stale 配置。

#### 3.5.2 `cc-ds.sh` 安全问题仍未修复（May26 已指出）

- 仍包含明文 API key：`ANTHROPIC_AUTH_TOKEN=sk-27ce3e9b...`。
- 仍使用 `--permission-mode auto`。
- 虽然 `.gitignore` 排除了该文件，但 key 仍暴露在本地磁盘上。

---

## 四、优先级修复建议

| 优先级 | 问题 | 建议修复 |
|--------|------|----------|
| **P0 (立即)** | `exp-D/md/launch.sh` 调用 exp-C 脚本 + 非法参数 | 1. 改为调用 `"$REPO/exps/exp-D/md/run_md.py"`；2. 移除 `--prmtop/--inpcrd/--outdir`，改用 exp-D `run_md.py` 接受的参数格式；3. `--gpu 0` 改为 `--gpu $gpu` |
| **P0 (立即)** | `build_complex.py` 化学错误的 LYS→ARG 突变 | 1. 停止使用原子重命名策略；2. 改用 PyMOL/Modeller 做侧链替换；3. 或使用 tleap `mutate` 命令在序列层面处理；4. 重新生成所有 exp-D 拓扑 |
| **P0 (立即)** | `build_complex.py` 生成的双 bond NZ（`bond NZ to N and C`） | 1. 移除 `bond complex.116.NZ LYA.1.C`；2. 仅保留一个正确的酰胺键；3. 重新运行所有 tleap 构建 |
| **P1 (本周)** | `common/scripts/run_md.py` barostat 在 Simulation 后添加 | 1. 将 `system.addForce(MonteCarloBarostat(...))` 移到 `Simulation(...)` 之前；2. 或添加 `simulation.reinitializeContext()` |
| **P1 (本周)** | 硬编码绝对路径 | 1. 统一改为 `pathlib.Path(__file__).resolve().parents[N]`；2. 在 `common/scripts/` 中引入 `paths.py` 统一路径管理；3. 所有 `run_md.py` 通过 `--repo` 参数接收根目录 |
| **P1 (本周)** | `addIonsRand Na+ 0 Cl- 0` 不加离子 | 1. 改为 `addIonsRand complex Na+ 0`（自动中和）；2. 或显式计算净电荷后添加对应数量的反离子 |
| **P1 (本周)** | exp-D 缺少文档 | 1. 编写 `README.md` 说明科学目标、体系、MD 方案；2. 编写 `exp-log.md` 记录构建过程和已知问题；3. 编写 `tasks.md` 跟踪 replica 进度 |
| **P2 (生产前)** | `aib_residue.xml` 缺少 Angle/Dihedral | 从 `aib.frcmod` / `aib_capped.frcmod` 提取 angle、proper dihedral、improper dihedral 参数，补充到 XML 模板中 |
| **P2 (生产前)** | `mmpbsa.py` 不完整 | 完成汇总表格的数据行输出；或从生产目录移除 |
| **P2 (生产前)** | 旧版脚本归档 | 将 `convert_to_aib8.py`、`fix_aib8_v2.py`、`build_diacid.py` 等旧版本移入 `archive/_versions/` |
| **P2 (生产前)** | 分析脚本统计方法升级 | 增加有效样本量校正、非高斯分布的中位数/IQR 报告、自相关时间估计 |
| **P3 (持续)** | `cc-ds.sh` API key 安全 | 1. 在 DeepSeek 控制台撤销该 key；2. 使用 `.env` 文件 + `python-dotenv` 注入；3. 考虑 `git filter-branch` 从历史中移除 |
| **P3 (持续)** | exp-C 假说修正 | 在 `exp-C/README.md` 和 `roadmap.md` 中正式修订"双点锚定"假说，将 linker-FA 逃逸作为新的科学发现而非失败记录 |

---

## 五、结语

这是一个**科学动机极强、执行力极高**的项目。在短短 5 天内，团队完成了从 exp-A 的 200 ns 生产到 exp-C 的 9 条 100 ns 验证，再到 exp-D 的 5 种 linker 变体构建与启动。exp-C 中 linker-FA 从 HSA 逃逸的 unexpected 发现，恰恰证明了计算模拟的价值——它揭示了仅靠游离 FA 实验无法预见的 linker 效应。

然而，**快速推进的代价是工程质量的显著恶化**。exp-D 的 launcher 调用错误脚本、build_complex.py 的化学错误突变、以及双 bond NZ 的化学荒谬性，都是"如果按原样运行就会得出错误数据"的 P0 级问题。这些 bug 的存在，使得 exp-D 当前产出的所有轨迹和拓扑文件都需要在修复后重新验证。

项目的 `best-practice.md` 第 35 条写着 "frame-0 验证：clash check + 关键距离 + 初始能量"。当前的 exp-D 体系恰恰需要这一验证：在启动生产 MD 前，必须确认 (1) Lys34→Arg 突变的几何是否合理，(2) LYA-linker 的酰胺键连接是否正确，(3) 初始能量是否收敛到合理范围。

**建议立即暂停 exp-D 的生产 MD**，修复上述 P0/P1 问题后，重新构建拓扑并做 frame-0 验证，再决定是否继续。否则，~7.5 μs 的 GPU 时间可能基于物理上错误的体系。

---

*Review 完成时间：2026-06-01*
*Reviewer：Kimi*
