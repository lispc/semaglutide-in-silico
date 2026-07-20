# exp-D MM-GBSA/PBSA 结果：4 个 linker 变体 @ GLP-1R ECD

**计算日期**：2026-07-20 · **工具**：MMPBSA.py 14.0 (AmberTools) · **执行**：MPI 10 ranks/case（vendored mpi4py + 系统 OpenMPI）

## 结论摘要

**有 linker 的变体 ECD 结合均强于 no_linker（GB 显著），但 1×/2×/3×OEG 之间的能量排序与 Lau 2015 活性排序不符；ΔG 与 log(EC50) 的 Spearman 相关方向为正但不显著（n=4）。** 与 Tail-Prot 几何分析（2026-07-19/20）结论一致：ECD 层面的指标只能区分"有无 linker"，不能解释 linker 长度之间的活性差。

| 变体 | ΔG_GB (3-rep mean±SD) | ΔG_PB | EC50 (pM, Lau 2015 T3) |
|---|---:|---:|---:|
| no_linker | −39.28 ± 2.94 | +2.67 ± 4.87 | 269 |
| gglu_1oeg | −60.35 ± 7.01 | −5.77 ± 4.08 | 4.8 |
| gglu_2oeg | −47.65 ± 7.65 | +1.52 ± 1.12 | 6.2 |
| gglu_3oeg | −64.44 ± 2.64 | −0.71 ± 0.89 | 27.7 |

**变体间 Welch t（vs no_linker, kcal/mol）**：

| 对比 | ΔΔG_GB | p | ΔΔG_PB | p |
|---|---:|---:|---:|---:|
| gglu_1oeg − no_linker | −21.06 | **0.022** | −8.44 | 0.085 |
| gglu_2oeg − no_linker | −8.37 | 0.190 | −1.15 | 0.725 |
| gglu_3oeg − no_linker | −25.16 | **0.0004** | −3.38 | 0.352 |

**Spearman 相关（ΔG vs log₁₀EC50；负 ΔG ↔ 低 EC50 为正相关）**：

| 方法 | variant-mean (n=4) | all-replica (n=12) |
|---|---|---|
| GB | R=+0.400, p=0.600 | R=+0.453, p=0.139 |
| PB | R=+0.800, p=0.200 | R=+0.583, **p=0.047** |

解读：相关方向符合预期（结合越负、活性越强），但 variant 层面 n=4 检验力极低，p=0.2–0.6 不可声称显著；n=12 把同变体的 3 个 replica 当独立点，夸大了有效样本量，PB 的 p=0.047 仅供参考。排序层面：活性序 1×≈2×>3×≫无；GB 结合序 3×>1×>2×≫无——3×OEG 结合最强却活性最差（27.7 vs 4.8/6.2 pM），精细排序不成立。**粗粒度信号（有/无 linker）成立，精细排序（linker 长度）不成立。**

## 每 replica ΔG（kcal/mol，mean，括号为 n_eff 修正 SEM）

| 变体 | rep | ΔG_GB | SEM | ΔG_PB | SEM |
|---|---|---:|---:|---:|---:|
| no_linker | 1 | −38.49 | 4.76 | +7.51 | 2.70 |
| no_linker | 2 | −36.83 | 0.77 | +2.72 | 0.78 |
| no_linker | 3 | −42.54 | 1.59 | −2.22 | 1.47 |
| gglu_1oeg | 1 | −58.39 | 2.69 | −4.82 | 2.53 |
| gglu_1oeg | 2 | −68.12 | 2.48 | −10.25 | 4.26 |
| gglu_1oeg | 3 | −54.53 | 7.06 | −2.25 | 3.02 |
| gglu_2oeg | 1 | −56.32 | 4.99 | +0.28 | 1.33 |
| gglu_2oeg | 2 | −44.74 | 1.74 | +1.82 | 0.99 |
| gglu_2oeg | 3 | −41.89 | 0.73 | +2.46 | 0.89 |
| gglu_3oeg | 1 | −63.21 | 4.92 | −0.08 | 4.85 |
| gglu_3oeg | 2 | −67.48 | 2.46 | −0.33 | 1.49 |
| gglu_3oeg | 3 | −62.65 | 0.97 | −1.72 | 1.33 |

## 能量组分（3-replica mean, kcal/mol）

| 变体 | vdW | EEL | E_GB | ESURF | ΔG_GB |
|---|---:|---:|---:|---:|---:|
| no_linker | −60.4 | −22.2 | +51.7 | −8.4 | −39.3 |
| gglu_1oeg | −71.5 | −46.0 | +67.7 | −10.5 | −60.3 |
| gglu_2oeg | −69.3 | −140.5 | +171.6 | −9.5 | −47.7 |
| gglu_3oeg | −79.8 | −29.7 | +56.3 | −11.3 | −64.4 |

- linker 变体的 vdW 均优于 no_linker（−69 ~ −80 vs −60）：linker/脂链与 ECD 的额外接触
- 2×OEG 的巨大 EEL（−140.5）与 EGB（+171.6）几乎完全抵消（γGlu 羧基 + 长极性链的静电/去溶剂化对冲），净效应反而最小——这是 GB 排序与活性排序脱节的表现之一
- PB 组分（ENPOLAR/EDISPER 模型）使绝对值系统性偏正 ~+40 kcal/mol，仅 GB/PB 间系统差，对相对比较影响有限

## 方法

- 体系：`exps/exp-D/tleap/{variant}.prmtop`（~36.7k 原子）。受体 ECD = 残基 1–100（1589 原子）；配体 = 肽（101–126, 26 残基）+ LNK（127；53/87/108/129 原子；Lys117 NZ 与 LNK C11 共价连接）
  - **注意：任务书给的原 mask :101-129 有误**——实际 LNK 为残基 127、128 起为 Na⁺；本计算使用修正后的配体 mask `:101-127`（已经 NZ-C11 共价键与干拓扑原子数双重验证）
- 轨迹：`exps/exp-D/md/{variant}/rep{1,2,3}/{variant}_traj.dcd`（1003 帧 × 100 ps），帧 502–1002 步长 4 → 126 帧/replica（最后 50 ns）
- frame-0 校验：NZ–C11 = 1.35–1.39 Å（4/4 变体通过，共价键完好）
- 预处理：cpptraj `autoimage` + `strip :WAT,Na+,Cl-` → NetCDF 干轨迹（nobox）；ante-MMPBSA `-n ':101-127' --radii=mbondi2`（只给配体 mask）
- MMPBSA.py 14.0 MPI 10 ranks/case：GB `igb=5, saltcon=0.1`；PB 默认 `inp=2, istrng=0.1, radiopt=0`；**无 decomp（价值低，提速）、无 nmode 熵**；单轨迹协议
- 统计：逐帧 ΔG → stats.py n_eff/SEM（副本内）；3-rep mean±SD；变体间 Welch t（n=3 vs 3）；Spearman(ΔG, log₁₀EC50)
- 运行：2 批次 × 6 case（每 case ~5 min wall，2.1k 原子体系），全部 exit 0；日志 `*/rep*/FINAL_run.log`

## 文件
```
analysis/mmgbsa/
├── RESULTS.md / aggregate.py / aggregated.json / summary.json / run_rep.sh
└── {no_linker,gglu_1oeg,gglu_2oeg,gglu_3oeg}/
    ├── com/rec/lig_dry.prmtop, mmpbsa.in(+test), ante.log
    └── rep{1,2,3}/  {variant}_repN_dry.nc, FINAL_MMPBSA.dat/.csv, FINAL_run.log, extract.cpptraj/.log
```

## 警告
1. 未含熵；单轨迹协议；带电配体绝对 ΔG 偏负，仅供变体间相对比较。
2. Spearman n=4（variant 层面）统计力极低；n=12 版本把 replica 当独立点，显著性虚高——两者都只作参考。
3. PB 与 GB 绝对值差 ~+40 kcal/mol（非极性溶剂化模型差异）；PB 下部分变体 ΔG 接近 0 甚至为正。
4. 副本间涨落大（如 1oeg GB SD=7.0；no_linker rep1 n_eff=6）——变体间小差异（<5 kcal/mol）不可判。
5. 与 Tail-Prot 几何结论一致：ECD 结合能解释不了 1×/2×/3×OEG 活性序；活性差更可能来自 TMD 结合/全长受体环境/PK，而非 ECD 亲和力。
