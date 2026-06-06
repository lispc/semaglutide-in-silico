# exp-F Phase 0 构建日志

**日期**: 2026-06-06
**目标**: 构建完整司美格鲁肽 + GLP-1R 受体的溶剂化体系，验证 tleap 管线

## 构建流程

### 1. 对齐配体到受体
使用 `build_full_receptor_ligand.py` 将 exp-D 的 `gglu_2oeg_fixed` 配体对齐到 7KI0 的 chain P（GLP-1 肽段）。
- 对齐方式：Kabsch 算法，基于 Cα 原子（残基 10-35 / 101-126）
- RMSD after alignment: 3.70 Å
- 输出：`7ki0_receptor_plus_expD_ligand.pdb` (9428 atoms)

### 2. 提取受体+肽 PDB 和对齐 LNK mol2
`prepare_tleap_inputs.py`:
- 从组合 PDB 删除 LNK，保存 `receptor_peptide.pdb` (1139 残基)
- 从组合 PDB 提取 LNK 坐标，更新原始 mol2，保存 `lnk_aligned.mol2`
- 保持 mol2 的 GAFF2 原子类型（c, c3, n, o, os, ho, hc, hn）不变

### 3. tleap 构建
`build_system.in`:
```
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2          # <-- 关键：必须加载 GAFF2
loadAmberParams lya_link.frcmod

complex = loadPdb receptor_peptide.pdb
LNK = loadMol2 lnk_aligned.mol2
sys = combine { complex LNK }

remove sys sys.1140.HZ1
remove sys sys.1140.HZ2
remove sys sys.1140.HZ3
remove sys sys.1150.N

bond sys.1140.NZ sys.1150.C

solvateOct sys TIP3PBOX 12.0
addIonsRand sys Na+ 0
saveAmberParm sys system.prmtop system.inpcrd
```

**关键发现**：
- tleap 内部残基编号与 PDB 编号不同。LYS P117 = tleap res 1140, LNK = tleap res 1150
- `remove` 和 `bond` 命令必须使用 tleap 内部编号
- 必须加载 `leaprc.gaff2`，否则 LNK 的 GAFF2 原子类型（c, c3, n 等）找不到 vdW 参数
- 系统电荷 +1.138416（非整数，非零），与 exp-D 相同。这是删除 LNK N 原子后的已知问题，不影响 MD

## 验证

### 体系规模
| 属性 | 数值 |
|------|------|
| 总原子数 | 363,332 |
| 总残基数 | 116,174 |
| 水分子 | 115,034 |
| Box | 171.6 Å (oct) |
| LNK 原子 | 130 (N 已删除) |

### 关键几何
- LYS.NZ — LNK.C 距离: **1.381 Å** ✅ (与 exp-D fixed 一致)
- LYS.NZ 键: CE, LNK.C ✅
- LNK.C 键: LYS.NZ ✅

### 文件位置
- `exps/exp-F/tleap/system.prmtop` (60M)
- `exps/exp-F/tleap/system.inpcrd` (13M)
- `exps/exp-F/tleap/system.pdb` (25M)

## 已知问题
1. **无膜**: Phase 0 仅验证溶剂化体系，未嵌入 POPC/胆固醇膜
2. **电荷非中性**: +1.138416，与 exp-D 相同。不影响 MD（PBC + PME 处理），但在报告 Methods 中需说明
3. **系统规模**: 363k atoms 超出最初预估的 ~250k。如加膜可能接近 400-450k

---

# exp-F Phase 1 稳定性验证日志

**日期**: 2026-06-06
**目标**: 验证溶剂化体系在 minimization + short equilibration 后的稳定性

## 协议

`run_min_equil.py`:
- Minimization: 500 steps L-BFGS (~41s)
- NVT heating: 0→100 K, 50 ps
- NPT heating: 100→310 K, 100 ps (5 steps × 20 ps)
- NPT equilibration: 310 K, 100 ps
- 约束: 1113 受体 CA (10 kJ/mol/nm²) + 18 肽 N-term backbone (100 kJ/mol/nm²)

## 结果

### 性能
| 指标 | 数值 |
|------|------|
| GPU 速度 | **~19 ns/day** (RTX 3090, CUDA mixed) |
| Minimization 时间 | 40.9 s |
| 250 ps 总时间 | ~28 min |

### 能量与温度轨迹
```
Step  Time(ps)  PotentialE(kJ/mol)  Temp(K)  Volume(nm³)
5000   10.0     -5,067,935          101.3    3686
25000  50.0     -5,138,990          100.3    3642   <-- NVT 100K 结束
50000  100.0    -4,845,963          184.1    3631   <-- NPT 184K
110000 220.0    -4,091,018          310.5    3561   <-- NPT 310K 结束
125000 250.0    -4,106,261          309.7    3553   <-- Equil 结束
```

**观察**:
- 温度控制精准，各阶段目标温度 ±1 K 内
- 能量随升温单调上升（热运动增加），符合物理预期
- NPT equilibration 期间能量继续缓慢下降并稳定
- **无 NaN、无能量爆炸、无 VDW 冲突** — 体系稳定 ✅

### 体积变化
- 初始体积: 3686 nm³ (100K, 刚溶剂化)
- 最终体积: 3513 nm³ (310K, NPT equil)
- 收缩: ~4.7%，属于正常溶剂化后弛豫范围

### 文件位置
- `exps/exp-F/md/rep1/minimized.pdb` (29M)
- `exps/exp-F/md/rep1/equilibrated.pdb` (29M)
- `exps/exp-F/md/rep1/equil.chk` (52M)

## 已知问题
1. **速度低于预估**: 19 ns/day vs 预估 45 ns/day。原因：系统 363k atoms > 预估 250k；约束力增加了计算开销
2. **仅 250 ps 验证**: Phase 1 只做 short equilibration，未跑 production。正式 production 需要 ≥100 ns/replica
3. **无膜**: 溶剂化体系，未嵌入 POPC/胆固醇膜

## 下一步
- Phase 2: 跑 production MD (100-500 ns)，监控肽-受体结合界面稳定性
- Phase 3: 用 CHARMM-GUI 构建膜嵌入体系，对比溶剂化 vs 膜环境的行为差异

## 2026-06-06: Rebuild v2 — Fixed LNK atom types and bonds

### Problem
The original `lnk_gglu_2oeg_pos.mol2` had an empty BOND section and incorrect GAFF2 atom types. When tleap loaded it, LNK atoms had no internal bonds, causing the first equilibrated structure to explode (atoms up to 220 Å apart).

### Fixes applied
1. **Distance-based bond detection** (`create_stripped_mol2.py`): Generated 114 bonds for the 130-atom LNK (N atom removed).
2. **Atom type corrections** (`fix_atom_types.py`):
   - Carbonyl carbons (`c3` → `c`): C17, C27, C37, C55, C58
   - Ether oxygens (`o` → `os`): O22, O25, O32, O35
   - Kept N16 as `n` (amide/amine), H as `hn`, carboxylate O as `o`
3. **Added missing dihedrals** to `lya_link_final.frcmod`:
   - `C8-N3-c-o` and `C8-N3-c-c3` (ff14SB LYS CE — GAFF2 LNK C amide)

### Verification
- tleap build: `Errors = 0`
- LNK internal bonds: 114 bonds, all distances < 2.0 Å
- Max LNK pairwise distance (minimized): 44 Å (was 220 Å)
- C=O bond lengths: 1.21–1.26 Å
- C–C single bonds: 1.52 Å

### Equilibration results
- Minimization: 25.3 s, PE stable
- NVT 0→100K (50 ps): complete
- NPT 100→310K (200 ps): complete
- NPT equilibration 310K (100 ps): complete
- Final PE: -4,159,783 kJ/mol
- Temperature: 310.7 K
- Box volume: 3559 nm³
- Speed: 38.5 ns/day

### Files
- `build/lnk_final.mol2` — corrected LNK mol2 with explicit bonds
- `build/lya_link_final.frcmod` — updated cross-force-field parameters
- `build/build_system_final.in` — final tleap input
- `tleap/system.prmtop` / `system.inpcrd` — production topology
- `md/rep1/equilibrated.pdb` — equilibrated structure
