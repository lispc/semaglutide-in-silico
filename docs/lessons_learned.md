# Semaglutide In-Silico 项目：经验教训与错误记录

> 本文档记录项目中的关键错误、误诊、及其根本原因，以防止重复发生。

---

## 错误 #1：Lys20–LNK 共价键"缺失"误诊（2026-06-06）

### 症状
- tleap 日志显示 `bond SYS.117.NZ SYS.129.C` 报错
- 后续 `bond SYS.145.NZ SYS.157.C` 成功
- **错误结论**：`system_ecd_v2.prmtop` 中 Lys20 与 LNK 之间无共价键

### 影响
- 基于"键缺失"假设，错误解释了 100ns MD 的高 RMSD（35–55 Å）为"结构不稳定"
- 报告 5.7–5.9 节全部基于错误前提
- 差点决定"重建体系并重跑 MD"

### 根本原因
1. **tleap 内部编号与最终 prmtop 编号的混淆**：`SYS.145.NZ` 和 `SYS.157.C` 是 tleap 在添加溶剂前的内部顺序编号，不是最终 prmtop 的原子编号
2. **未验证即下结论**：看到 tleap 日志中第一行 bond 报错，就假设键未建立，没有用 `parmed` 或 MDAnalysis 检查最终 prmtop
3. **bond 命令的成功执行被忽略**：日志中明确显示第二行 bond 成功，但未意识到它对应的是同一对原子

### 验证方法（正确做法）
```python
import parmed as pmd

prm = pmd.load_file('system_ecd_v2.prmtop')
# 检查特定原子对是否成键
for bond in prm.bonds:
    if (bond.atom1.name == 'NZ' and bond.atom2.name == 'C') or \
       (bond.atom1.name == 'C' and bond.atom2.name == 'NZ'):
        print(f"Bond: {bond.atom1.residue.name}{bond.atom1.residue.number}.{bond.atom1.name} -- "
              f"{bond.atom2.residue.name}{bond.atom2.residue.number}.{bond.atom2.name}")
        print(f"  k = {bond.type.k:.1f} kcal/mol/Å², req = {bond.type.req:.3f} Å")

# MDAnalysis 二次验证
import MDAnalysis as mda
u = mda.Universe('system_ecd_v2.prmtop')
print("Bonds in topology:", u.bonds is not None)
print("Number of bonds:", len(u.bonds) if u.bonds else 0)
```

### 教训
1. **永远验证最终文件，不要只看中间日志**
2. **tleap 的内部编号 ≠ 最终 prmtop 编号**，bond 命令的成功与否应以最终拓扑为准
3. **多个独立工具交叉验证**：parmed + MDAnalysis + 可视化（VMD）三重确认

---

## 错误 #2：superposition=False RMSD 误导诊断（2026-06-06）

### 症状
- ECD v2 100ns MD 的 RMSD 高达 35–55 Å
- **错误结论**："结构不稳定，仍在漂移"

### 影响
- 与错误 #1 叠加，形成了"键缺失导致结构不稳定"的错误叙事
- 掩盖了真实的物理图景：ECD/肽段/LNK 内部结构稳定，通过柔性 linker 做域间运动

### 根本原因
- **未区分 RMSD 的两种含义**：
  - `superposition=True`：先最佳叠合再计算 RMSD → 衡量**内部构象变化**
  - `superposition=False`：直接计算坐标差 → 包含**平移/旋转/构象变化**
- 对于含柔性 linker 的体系，superposition=False 会将正常的域间运动（几十 Å）误判为"结构破坏"

### 正确分析方法
```python
from MDAnalysis.analysis import rms, align

# 1. 内部稳定性：superposition=True
# 适用于评估蛋白质/肽段是否保持折叠状态
rmsd_internal = rms.rmsd(sel.positions, ref_sel.positions, superposition=True)

# 2. 域间运动：质心距离
# 适用于评估 linker 柔性和域间相对位置
def com_distance(u, sel1, sel2):
    return np.linalg.norm(sel1.center_of_mass() - sel2.center_of_mass())

# 3. Linker 柔性：半径回转 Rg
rg = sel.radius_of_gyration()
```

### 教训
1. **柔性 linker 体系必须同时报告**：
   - 内部 RMSD（superposition=True）→ 结构稳定性
   - COM 距离 → 域间运动幅度
   - Rg → linker 伸展程度
2. **RMSD 数值大 ≠ 结构不稳定**，必须结合 superposition 参数和物理图景解读
3. **分析前先思考**：该体系的预期行为是什么？对于 HSA-linker-ECD 体系，大尺度域间运动是设计意图

---

## 错误 #3：分析脚本性能问题（2026-06-06）

### 症状
- `analyze_ecd_v2_100ns.py` 运行超时（>30 分钟）
- 4 个独立的 `AlignTraj` 调用，每个都要重新加载 100ns DCD（~4GB）

### 根本原因
- **重复加载大文件**：每个 RMSF 分析都创建新的 Universe 并重新加载 DCD
- **in_memory=False**：每次对齐都写入临时文件，I/O 开销巨大
- **未复用已对齐的轨迹**：ECD、Peptide、HSA 可以用同一个全局对齐的轨迹计算 RMSF

### 优化方案
```python
# 优化前（慢）：4 次独立对齐
for sel in ['resid 1-100', 'resid 101-128', 'resid 130-711', 'resid 129']:
    u = mda.Universe(PRMTOP, DCD)
    align.AlignTraj(u, ref, select=sel, in_memory=False).run()
    RMSF(...).run()

# 优化后（快）：1 次全局对齐，复用
u = mda.Universe(PRMTOP, DCD)
align.AlignTraj(u, ref, select='protein and name CA', in_memory=True).run(step=STEP)
# 所有 RMSF 都在已对齐的 u 上计算
ecd_rmsf = RMSF(u.select_atoms('name CA and resid 1-100')).run(step=STEP)
pep_rmsf = RMSF(u.select_atoms('name CA and resid 101-128')).run(step=STEP)
hsa_rmsf = RMSF(u.select_atoms('name CA and resid 130-711')).run(step=STEP)
```

### 教训
1. **I/O 是大 MD 分析的瓶颈**：避免重复加载轨迹文件
2. **`in_memory=True`**：当内存足够时（4GB DCD → ~8GB 内存），使用内存对齐避免临时文件
3. ~~**全局对齐优于局部对齐**：对多组分体系，用共同的骨架对齐一次，再分别计算各组分的 RMSF~~
   > ⚠️ **SUPERSEDED（2026-06-09）**：本条被后续实践推翻。多组分复合物中全局对齐会把域间运动计入各组分的 RMSF，正确做法是**每个组分独立对齐**（component-wise alignment）。见 `best-practice-v2.md` #18 及 exp-F ECD v2 的修正分析。上文"根本原因"中"可以用同一个全局对齐的轨迹计算 RMSF"同样作废；本条保留仅作为性能优化（避免重复加载、`in_memory=True`、`step` 跳帧）的参考。
4. **`step` 参数**：分析时跳过帧（如 step=5）可大幅加速，不影响统计可靠性

---

## 方法论总结

### 构建阶段
| 检查项 | 工具 | 何时做 |
|--------|------|--------|
| 共价键是否存在 | `parmed` / MDAnalysis | 构建后立即 |
| 原子编号对应关系 | `parmed` 打印 residue/atom 列表 | 困惑时 |
| 结构可视化 | VMD / PyMOL | 运行 MD 前 |
| 体系电荷中性 | `parmed` 统计 | 构建后 |

### 分析阶段
| 体系类型 | 必报指标 | 禁止做法 |
|----------|---------|---------|
| 柔性 linker 体系 | 内部 RMSD + COM 距离 + Rg | 仅报 superposition=False RMSD |
| 膜蛋白体系 | 膜法向对齐 + z-position COM 校正 | 绝对坐标 |
| 结合态复合物 | 界面距离 + 氢键 + MM-GBSA | 仅报整体 RMSD |

### 诊断流程
1. **数据异常？先怀疑工具/方法，再怀疑物理**
2. **大 RMSD？先问 superposition=True 还是 False**
3. **键是否建立？用 parmed + MDAnalysis 双重验证，不要只看 tleap 日志**
4. **性能问题？检查 I/O 和内存使用，优先复用加载的轨迹**

---

## 错误 #4：膜体系 LNK mol2 含垃圾原子导致 MD 能量灾难（2026-06-06）

### 症状
- 膜体系生产 MD（68 ns）中 PE 出现大量跳变：72% 的帧间差 > 1000 kJ/mol
- LNK 残基（resid 1140）含 47 个原子而非预期的 43 个
- VMD 可视化显示 LNK 链末端有 4 个孤立原子（C58/O59/N60/C61）在 ~40Å 外飘浮
- `parmed` 检查：`lnk_noh_fixed.mol2` 中无 C55-C58 键，C58 与主链完全断开

### 影响
- **68 ns MD 数据作废**：LNK 内部几何完全错误，脂肪酸链的物理行为不可信
- **受体/肽段/膜的分析仍可用**：TM6 位移、肽-受体相互作用、膜性质等不依赖 LNK 正确性
- **项目进度延迟**：需重建膜体系（重新 minimization + equilibration + production）

### 根本原因
1. **mol2 文件未经验证直接用于膜体系构建**：`lnk_noh_fixed.mol2`（以及 `lnk_dum.mol2`、`lnk_noh.mol2`）均含 4 个额外原子
2. **构建流程未检查原子连通性**：tleap 加溶剂、加膜时未报错，因为额外原子只是"孤立"而非化学上不可能
3. **Mol2 原子顺序问题**：mol2 中原子名与 PDB 中的 `ATOM`/`HETATM` 命名不同，导致很难直观对比

### 验证方法（正确做法）
```python
import parmed as pmd

# 1. 检查原子数
prm = pmd.load_file('lnk_noh_fixed.mol2')
print(f"Atoms: {len(prm.atoms)}")  # 应为 43，实际是 47

# 2. 检查末端连通性
for bond in prm.bonds:
    if 'C55' in [bond.atom1.name, bond.atom2.name]:
        print(f"C55 bonded to: {bond.atom1.name if bond.atom2.name == 'C55' else bond.atom2.name}")
# 正确：C55 应连 C54 和 C56（或 O56，取决于命名）
# 实际：C55 只连 C54，C58 未与 C55 成键

# 3. 可视化检查
# VMD: 加载 system.pdb，用 "resid 1140" 选择 LNK，旋转查看链末端
# 若有原子远离主链 → 立即停跑
```

### 修复步骤
1. **定位问题**：比较 LNK mol2 与正确结构（文献 SMILES / 原始构建脚本）
2. **移除垃圾原子**：用文本编辑器删除 mol2 中 4 个孤立原子的记录
3. **验证修复后**：
   - 原子数 = 43
   - 键数 = 42（或 43，取决于末端羧基）
   - VMD 可视化确认链连续
4. **重建整个膜体系**：从 `minimization` 开始，不能用旧的 `inpcrd`

### 教训
1. **mol2 是小分子拓扑的薄弱环节**：mol2 格式宽松，无标准验证工具，容易混入错误原子
2. **任何新的小分子 mol2 必须通过三重检查**：
   - 原子数 vs 预期 SMILES
   - 末端连通性（parmed 键列表）
   - 可视化（VMD/PyMOL）
3. **长链脂肪酸 mol2 特别脆弱**：C18 二酸链长，原子编号容易出错，末端几个原子最容易遗漏或重复
4. **能量异常是拓扑错误的第一信号**：PE 跳变 > 1000 kJ/mol 不应被忽视，应立即停跑检查
5. **构建脚本应加入自动化验证**：在 `tleap` 后、MD 前，检查所有非标准残基的原子数和连通性

---



---

## 2026-07-25 — exp-G：单轨迹结论在多稳态体系里不可信（n=3 是下限）

### 现象
exp-G（HSA+司美格鲁肽+ECD 三元体系）pilot 20 ns 显示"预锚定 FA3 ~2 ns 脱锚、尾羧基保持 HSA 表面接触"，
看似桥接构型的定论。production 3×100 ns 驻留分区分析却发现**三个 replica 三种命运**：
rep1 全程 FA3↔表面动态交换 346 次不逃离（FA3 占 55.6%），rep2 3.9 ns 脱锚转表面，rep3 2.5 ns 脱锚渐溶剂化。
pilot 的 20 ns 只代表 rep2/3 型行为——若只有 pilot，会完全错过"FA3 快速交换驻留"这一主要模式。

### 教训
1. **多稳态/缓慢交换体系，单轨迹任何时长的结论都可能是少数派样本**：驻留分区、脱锚时间这类指标必须多 replica 分别统计再合并，合并 SD 巨大本身就是信号
2. **pilot 只能用于烟测协议/速度，不能用于科学结论**：pilot 发现要写进文档时必须标注"单轨迹，待多 replica 验证"
3. **构建几何 ≠ 平衡几何**：exp-G frame-0 的近端 O38-ARG408 盐桥（放置时扫到 3 Å）在 3 个 replica 中成桥率 0%——frame-0 校验通过只说明无 clash、键完整，不证明相互作用在平衡态存在
4. **度量口径要先定义再报数**：pilot 条目里"肽-ECD 1.34 Å"度量不明无法复现（production 用最小重原子距离口径 frame-0=2.92 Å）——报界面距离时必须写清楚是"最小重原子距离"还是其他口径

*最后更新：2026-07-25*
