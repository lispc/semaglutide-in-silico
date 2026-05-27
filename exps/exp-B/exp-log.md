# exp-B 实验日志

> 只追加，不删除。每次记录含日期、时间、操作内容和结果。

---

## 2026-05-27 — 实验目录初始化

- 创建 `exps/exp-B/` 目录结构
- 编写 README.md（实验设计）、tasks.md（任务清单）、exp-log.md（本文件）
- **状态**：设计阶段，等待用户确认后启动 Phase 0

---

## 2026-05-27 — Phase 0: 结构准备（重大发现）

### 关键结构分析

下载并分析了三个结构：

| PDB | 类型 | Chain A (ECD) | Chain B/P (Peptide) | ECD-Peptide 最近距离 |
|-----|------|---------------|---------------------|---------------------|
| 3IOL | X-ray 2.10 Å | 残基 29–128 (WT) | 残基 10–37 (WT, K34) | **36.1 Å** |
| 4ZGM | X-ray 2.70 Å | 残基 29–128 (WT) | 残基 10–37 (R34, sema骨架) | **24.9 Å** |
| 7KI0 | cryo-EM 3.30 Å | GLP-1R 全长 (残基 29–423) | 残基 7–36 (sema, R34) | **38.6 Å** |

**意外发现**：所有已发表结构中，ECD 与肽链在 PDB 坐标中均不直接接触（最近原子距离 25–39 Å）。这在晶体结构中可能是晶格堆积造成的，但在 cryo-EM 结构中也不接触则暗示活性态 GLP-1R 中 ECD 可能远离 TMD 正位点。

此外，文献报道的 Arg34–Glu27 水介导氢键中的 Glu27 在所有结构中均缺失（REMARK 465: 残基 24–28 disordered）。

### 拓扑构建

- 使用 tleap (ff14SB + GAFF2 + TIP3P) 成功构建 B1 (K34) 和 B2 (R34) 拓扑
- Key fix: disulfide cysteines 需重命名为 CYX（非 CYS），否则缺少 S-S 扭转参数导致 tleap 失败
- 系统大小：~36k atoms, ~11.5k 水分子，盒子 ~77×67×87 Å
- 产出：`b1.prmtop` (6.1 MB), `b2.prmtop` (6.1 MB)

### 肽起点与限制

- 晶体结构只解析了肽残基 10–37（残基 7–9 包括 Aib8 在所有结构中 disordered）
- 当前拓扑只包含残基 10–37（使用 NGLY/CGLY 末端），Aib8 不在模型内
- 不影响 K34R 问题（位置 34 在 C 端，远离缺失的 N 端）

### 搜索与决策

- ddgr 搜索了 6GB1 (GLP-1R ECD + Peptide 11)、5VEW/5VEX (TMD-only)、6X18 等替代结构
- 发现 GLP-1R ECD 晶体结构中 ECD-肽分离是**系统性特征**（所有检查的结构 min dist = 25–39 Å）
- 6GB1 (space group P43 21 2, 2.73 Å) 同样 36.6 Å — 不同空间群、不同肽，相同现象
- 确认 BIOMT assembly = identity matrix，无结晶学对称操作可修复
- **决策**：使用 PCA 几何对齐 + Rosetta FlexPepDocking refinement

### PCA 几何定位

- 计算肽螺旋 (残基 13–33) 的 PCA 主轴，对齐到 ECD 疏水沟槽 (残基 32/39/69/89/90/91/123/124/125/127/128) 的 PCA 主轴
- 对齐后 min ECD-peptide distance: 3.8 Å（从 24.9 Å 显著改善）
- 严重 clash: 34 对 <1.5 Å（从 338 对显著改善）
- 已启动 Rosetta FlexPepDocking refinement（15 decoys, pep_refine mode）

### Rosetta FlexPepDocking 失败

- 10 decoys 生成后，最优 decoy (0010, total_score=-190.6) 的 ECD-peptide min dist 仍为 21.9 Å
- FlexPepDocking 仅优化了肽和 ECD 的内部几何，未能克服 25 Å 的初始分离
- 符合 best-practice §13 的预期：FlexPepDock 对长肽 (>10 残基) 不可靠

### 决策：暂停 exp-B，转向 exp-C

- **2026-05-27**：用户决策跳过 exp-B，优先推进 exp-C (HSA 脂链结合)
- exp-B 的 K34R 问题可通过替代方案后续解决：
  - AF3/Boltz-1 预测 ECD-肽复合物
  - 或使用全长 GLP-1R cryo-EM 结构（需处理膜环境）
  - 或简化为 K34 vs R34 肽单体比较
- B1/B2 tleap 拓扑文件已保存，可随时恢复

---
*维护者：Claude Code*
*最后更新：2026-05-27*
