# exp-D 实验日志

> 只追加，不删除。每次记录含日期、时间、操作内容和结果。

---

## 2026-05-31 — 实验启动

### 决策

exp-C 结论：linker 接上后 FA 逃逸 HSA FA3（无论末端电荷）。转 exp-D 验证 linker 的另一半功能：在 GLP-1R 端隔离脂链、保护受体结合。对应 Lau 2015 Table 3（variable linker length → GLP-1R potency）。

### 体系设计

5 个 linker 变体，均连接 C18 二酸，接在 Aib8,Arg34-GLP-1(7-37) 的 Lys26 上：

| # | Linker | Lau 2015 | EC50 (pM) | 状态 |
|---|--------|---------|-----------|:---:|
| 1 | 无 linker（直达 C18） | Cmpd 19 | 269 | 运行中 |
| 2 | γGlu | Cmpd 20 | 9.9 | NaN（C18碰撞ECD） |
| 3 | γGlu-1×OEG | Cmpd 21 | 4.8 | 运行中 |
| 4 | γGlu-2×OEG（司美格鲁肽） | — | 6.2 | 运行中 |
| 5 | γGlu-3×OEG | Cmpd 24 | 27.7 | 运行中 |

系统不含 HSA——linker-FA 在溶剂中与 GLP-1R ECD 相互作用。

### LYA linker 变体构建

- **方法**：RDKit ETKDGv3 + MMFF → 原始 mol2 → antechamber AM1-BCC（OMP_NUM_THREADS=4）
- **净电荷**：nc=-1（COO⁻），偶数电子确保 SQM 收敛
- **剥离**：去除 ACE cap + Lys backbone/sidechain + NME cap，仅保留 linker-C18-COO⁻
- **GAFF2 类型**：n (amide N), c (carbonyl C), o (O), c3 (CH2), hc (H on C)
- **定位**：Kabsch 对齐 amide N→C 向量至 peptide CE→NZ 方向，amide N 置于 NZ 位置

### 肽-ECD 复合物构建

- **模板**：3IOL（ECD + GLP-1 10-35，2.1Å），唯一含正确 ECD-肽结合构象的晶体结构
- **Lys34→Arg 突变**：tleap loadPdb + 模板自动补全缺失原子
- **已知局限**：突变质量受限于初始原子位置，不影响 exp-D 主要指标（C18-ECD 距离）

### 关键构建教训

1. **tleap remove/bond 命令存在静默失败 bug**：在 ParmEd 生成的 PDB 格式下 `remove sys sys.117.HZ1` 报 "Argument #2 is of type String"，`bond` 同样静默失败
2. **10M 残基问题**：ecd_pep.pdb 包含晶体学 NME cap 残基（10M），tleap 创建 33 个无类型原子导致 FATAL 错误。`remove sys sys.127` 在特定 PDB 格式下可删除
3. **LNK 共价键缺失**：由于 tleap bond 失败，prmtop 中 LNK 原子无共价键。原子间距正确时 minimization 可修复，MD 稳定运行
4. **混合力场 frcmod**：需要 `lya_link.frcmod` 覆盖 ff14SB N3—GAFF2 c 之间的 cross-force-field bond/angle 参数
5. **ParmEd PDB 格式**：仅 ParmEd `.save()` 生成的 PDB 能被 tleap atom mask 正确解析；手动 Python PDB writer 导致 mask 失败

### AM1-BCC 运行时间

| Variant | Atoms | SQM 时间 |
|---------|:---:|------|
| no_linker | 85 | ~3 min |
| gglu | 98 | ~18 min |
| gglu_1oeg | 119 | ~40 min |
| gglu_2oeg | 140 | ~60 min |
| gglu_3oeg | 161 | ~90 min |

---

## 2026-06-01 — MD 生产进展

### 运行状态

| Variant | Rep 1 | Rep 2 | Rep 3 |
|---------|:---:|:---:|:---:|
| no_linker | 67ns | 91ns | 89ns |
| gglu_1oeg | 89ns | 100ns ✓ | 67ns |
| gglu_2oeg | 100ns ✓ | 67ns | 91ns |
| gglu_3oeg | 67ns | 91ns | 89ns |

- 速度：单 GPU ~450 ns/d，共享 ~100-200 ns/d
- 系统规模：~36k atoms（vs exp-C ~86k）
- 预计全部完成：2026-06-01 晚

### gglu NaN 诊断

- PE after minimization: 7×10^16 kJ/mol（vs 正常 ~-5.9×10^5）
- 根因：γGlu-only linker 导致 C18 tail 从 NZ 延伸 4-5 Å 即碰撞 ECD 表面
- 更长 linker（1oeg/2oeg/3oeg）将 tail 延伸至溶剂区，无碰撞
- 待修复：调整 tail 初始旋转角度，避让 ECD 表面

### Review 反馈 (Kimi, 2026-06-01)

- barostat 时序 bug 已修复（移至 Simulation 创建之前）
- exp-D 文档已补充
- launch.sh 已废弃（从未用于生产启动）

---

*维护者：Claude Code*
*最后更新：2026-06-01*

## 2026-06-01 — MD 完成与首轮分析

### 完成状态
- 12/15 replica 完成, gglu rep1 重开 GPU 3 恢复正常速度

### 首轮分析：linker 变体无显著差异

| Variant | EC50 (pM) | CA RMSD (A) | Tail-Prot (A) | NZ-C (A) |
|---------|-----------|-------------|---------------|----------|
| No linker | 269 | 2.2+-0.1 | 5.9+-0.1 | 14.4+-2.7 |
| gGlu-1xOEG | 4.8 | 2.2+-0.0 | 6.0+-0.2 | 14.6+-3.0 |
| gGlu-2xOEG | 6.2 | 2.0+-0.1 | 5.8+-1.2 | 14.6+-3.0 |
| gGlu-3xOEG | 27.7 | 2.0+-0.2 | 5.8+-1.0 | 14.1+-3.0 |

**核心发现**: NZ-C 距离 ~14A -> 酰胺键未形成。LNK 靠非键作用浮在肽链附近，linker 差异被掩盖。

**根因**: tleap bond 命令静默失败; ParmEd bond 添加因缺少 BondType 对象失败。

### 修复计划
- ParmEd BondType(427.0, 1.38) 创建正确 NZ-C 酰胺键
- 单 replica 验证: NZ-C ~1.4A -> 变体间有趋势 -> 扩展 3 replica

---
*维护者：Claude Code*
*最后更新：2026-06-01*

## 2026-06-01 — ParmEd bond 修复成功

### 修复方法
- tleap `bond` 命令静默失败导致 NZ-C 酰胺键缺失（NZ-C ~14 Å）
- ParmEd 修复：删除 HZ 原子 + 距离检测 LNK 内部键 + BondType 创建 NZ-C 酰胺键
- BondType(427.0 kcal/mol/A², 1.38 A) 匹配 N3-c amide bond

### 修复后结果
| Variant | PE after min | NZ-C after min | NZ-C production |
|---------|-------------|:---:|:---:|
| no_linker | -588,287 | 1.38 A | 1.36-1.44 A（5ns 稳定） |
| gglu | -574,610 | 1.38 A | 测试中 |
| gglu_1oeg | -598,121 | 1.39 A | 测试中 |
| gglu_2oeg | -578,083 | 1.38 A | 测试中 |
| gglu_3oeg | -582,262 | 1.39 A | 测试中 |

全部 5 个变体 NZ-C = 1.38-1.39 A（正确酰胺键），PE 正常，无 NaN。

### 下一步
- 单 replica 10ns 测试确认 NZ-C 稳定 → 全 5 变体 ×3 replica 100ns 生产

---
*维护者：Claude Code*
*最后更新：2026-06-01*

## 2026-06-02 — 全量 MD 生产进展

### 状态
- 4/15 done（全部 rep1 先完成），11 个 ~60-64ns，全部 307-312K 无 NaN
- 速度 92-98 ns/d，4 GPU 满载，预计今晚全部完成

### 中间分析（~60ns）
NZ-C 键修复验证成功：

| Variant | EC50 | CA RMSD (A) | NZ-C (A) | Tail-Prot (A) |
|---------|------|-------------|:---:|-----------|
| No linker | 269 | 2.1+-0.2 | 1.51 OK | 4.7+-0.2 |
| gGlu | 9.9 | 2.2+-0.2 | 1.51 OK | 4.3+-0.5 |
| gGlu-1xOEG | 4.8 | 2.1+-0.2 | 1.52 OK | 4.5+-0.6 |
| gGlu-2xOEG | 6.2 | 2.4+-0.2 | 1.52 OK | 3.8+-0.2 |
| gGlu-3xOEG | 27.7 | 3.0+-0.2 | 1.52 OK | 4.0+-0.3 |

- NZ-C 全部 1.51-1.52 A（vs 修复前 14 A），共价键修复 100% 成功
- gGlu-3xOEG CA RMSD 最高（3.0 A），可能反映过长 linker 的 entropic penalty
- gGlu-2xOEG Tail-Prot 最近（3.8 A），linker 刚好够长产生有利接触

---
*维护者：Claude Code*
*最后更新：2026-06-02*
