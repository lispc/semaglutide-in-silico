# exp-G: HSA–受体竞争三元复合物（预锚定 FA3）

> 对应核心问题：司美格鲁肽 99% 结合 HSA，受体激活由游离态完成；2×OEG 的真正选择依据是
> BR ratio（HSA 亲和力 vs 受体亲和力的平衡）。exp-D 只有 ECD 端（无 HSA），exp-C 只有
> HSA 端（无受体），exp-F 把三者放进一盒但 HSA 游离放置 111 Å、100 ns 从未接近。
> exp-G 的前提：**脂链尾预锚定到 HSA FA3**，直接观测两端的竞争。

## 体系组成

| 组分 | 来源 | 说明 |
|------|------|------|
| HSA | `exp-C/tleap/hsa_no_myr.pdb` | 582 残基（晶体编号 3–584，链 A），去肉豆蔻酸 |
| ECD + 肽 | `exp-D/tleap/ecd_pep_nocap.pdb` | ECD 1–100（链 R）+ GLP-1(10-35),Aib8,Arg34 26 残基（链 P），3IOL 结合构象 |
| LNK | `exp-D/tleap/lnk_gglu_2oeg_clean.mol2` | γGlu-2×OEG-C18 二酸（108 原子），GAFF2 类型+完整键表（exp-D 重建版） |
| 水/离子 | TIP3P，缓冲 12 Å（orthorhombic `solvateBox`，体系细长） | Na+/Cl- 中和 + 0.15 M NaCl |

原子数预算：溶质 ~12k（HSA 9.2k + ECD 1.6k + 肽 0.5k + LNK 0.1k）+ 水 ~55–65k → **总计 ~180–210k**。
tleap 编号（连续）：HSA 1–582，ECD 583–682，肽 683–708（Lys26 = 肽第 17 残基 = **699**），LNK = **709**。

## 预锚定策略（核心设计）

脂链按 exp-C MM-GBSA 验证的 FA3 几何放置（1E7G MYR 1003 晶体位姿）：

1. **远端羧基**（LNK C55/O56/O57，游离 COO⁻）置于 MYR C1 位，羧基平面按晶体 O1/O2 搭建
   → 与 ARG346+ARG483 形成双盐桥（frame-0 要求 O→胍基 ≤3.5 Å）。
2. **C18 链 C54–C45** 沿晶体 MYR C2–C11 穿越口袋；**C44–C39** 贪心二面角扫描延伸出口袋。
3. **近端羧基**（酰胺羰基 C37=O38，C18 二酸接 ADO 的一端）朝向 ARG408：
   放置时扫描使 O38 距 ARG408 胍基 ~3 Å（frame-0 ≤5 Å）。
4. **ADO2/ADO1/γGlu**（N29…C13）以贪心扫描长出口袋、伸向溶剂（clearance + 外向引导 + 紧致约束）。
5. **肽-ECD 复合体**作为刚体（3IOL 结合构象不动）放置：Lys26 NZ 落在 C11 酰胺 N 位
   （NZ-C11 = 1.38 Å，CE-NZ-C11 共线，与 exp-D 的 NZ-C11 构建一致）；绕 CE-NZ 轴扫描旋转，
   使 ECD/肽与 HSA 零 clash（<2.2 Å 重原子对 = 0）且 ECD 尽量朝外。
   初始 ECD-HSA 间距由 linker 自然伸展跨度决定（C37→NZ 共 22 个重原子，伸展 ~27 Å），
   允许大间距——**不做人为拉近**。

## 初始几何松弛（约束分级释放）

pilot/生产统一协议（`md/run_pilot.py`）：

| 阶段 | 时长 | HSA+ECD CA 约束 k | 锚点 O（O56/O57/O38）k |
|------|------|:---:|:---:|
| minimize | 10,000 步 | 5.0 | 5.0 |
| NVT 0→100 K | 50 ps | 5.0 | 5.0 |
| NPT 100→310 K | 100 ps | 5.0 | 5.0 |
| NPT eq | 200 ps | 5.0→1.0（eq 结束切换） | 5.0→0.0 |
| production | 310 K NPT | **1.0** | **0.0（完全释放，否则无法观测竞争）** |

肽与 LNK 全程不加约束。k 用 OpenMM global parameter 实现，生产中可切换。

## 进展（2026-07-24）

- **Pilot 20 ns ✅**：预锚定双盐桥在生产开始后 ~2 ns 松脱（linker 使 FA3 解离加速 ~100×），
  但尾羧基全程保持 HSA 表面接触（2.53 ± 0.14 Å，后 25%），肽-ECD 界面 1.34 ± 0.03 Å 稳定——
  **桥接构型是稳态**，γGlu-2×OEG 让分子同时挂住受体与白蛋白表面。详见 `exp-log.md` 2026-07-21 条目
- **Production**：rep1/rep2 各 100 ns ✅（seeds 101/202），rep3（seed 303）运行中（GPU 3 共享，~66 ns/d）
- **待办**：3 replica 齐后做驻留分区分析（FA3 口袋 <4 Å / HSA 表面 ~2.5 Å / 水相 >10 Å 帧占比 +
  近端羧基-ARG408 成桥率 + 肽-ECD 稳定性），结果写 `analysis/RESULTS.md`
- 已知小瑕疵：pilot 因 0-based 编号多约束了肽 N 端 1 个 CA，无实质影响；`run_pilot.py` 已修正（`<=681`），生产版生效

## 观测指标（分析脚本后续补）

- **锚定稳定性**：远端 COO⁻ O → ARG346/483 胍基最小距离；近端 O38 → ARG408 距离（逐帧）
- **肽-ECD 界面**：肽(683–708) 与 ECD(583–682) 最小重原子距离 + 界面接触数
- **LNK 形态**：重原子 Rg、C11→C55 端到端距离（跨度）
- **HSA–ECD 接近**：两蛋白 CA 最小距离、接触数
- **竞争事件判据**：远端锚断裂 = O56/57→Arg 距离 >6 Å 持续 >1 ns；肽端脱钩 = 肽-ECD 界面最小距离 >6 Å 持续 >1 ns；记录先发生的一端及时间

## Frame-0 验证判据（构建后必过）

1. 远端羧基 O → ARG346/ARG483 胍基 ≤3.5 Å
2. 近端酰胺 O38 → ARG408 ≤5 Å
3. 肽-ECD 界面最小距离 ≤4 Å
4. HSA↔ECD/肽/LNK 重原子 clash（<2.2 Å）= 0；LNK 自接触（1-4+）≥2.1 Å
5. NZ-C11 键存在（parmed 验证 BondType 427/1.38），0 个未参数化项
6. 总原子数/盒尺寸/离子数记录；净电荷 ~0
7. OpenMM minimize 500 步 PE ≈ -3e6 kJ/mol 量级（~200k 原子体系），无 NaN

## MD 协议（与 exp-D 对齐）

ff14SB + GAFF2 + TIP3P；PME 1.0 nm；HBonds 约束；2 fs；Langevin 310 K, 1 ps⁻¹；
MonteCarloBarostat 1 bar（**Simulation() 之前加入**）；DCD 每 50,000 步；log 每 10,000 步。

## 文件结构

```
exp-G/
├── README.md / exp-log.md
├── build/
│   ├── place_lnk_fa3.py      # LNK 尾部锚定 FA3 + 链生长 + NZ 目标计算
│   ├── assemble_complex.py   # 组合 PDB（parmed 写盘，避免 mask 解析坑）
│   ├── validate_frame0.py    # 上述全部 frame-0 判据
│   └── preflight_min.py      # OpenMM 500 步 PE 检查
├── tleap/
│   ├── lnk_2oeg_fa3.mol2     # 锚定几何 LNK（坐标新建，类型/电荷/键沿用 exp-D 重建版）
│   ├── complex_ecd_pep_hsa.pdb
│   ├── build_complex.in (+ _ions 两遍：先溶水后按水数加 0.15 M 盐)
│   └── complex.prmtop/.inpcrd
└── md/
    ├── run_pilot.py          # 约束分级释放版 run 脚本
    └── pilot/                # GPU3 20 ns 烟测输出
```

## 已知风险与对策

- **linker 跨度不足**（C37 锚定位姿与 NZ 目标距离 >30 Å，链被拉直）：不强行连接；
  若贪心生长末端偏离 >8 Å，重做锚定段二面角（记录于 exp-log）。
- **200k 原子 pilot 在共享 GPU 上慢**：烟测只要 20 ns 验证稳定性即可，不等跑完。
- **锚点释放后瞬时脱钩**：若 pilot 中远端 COO⁻ 在释放后 <5 ns 内 >6 Å，说明初始
  盐桥几何或 H  placement 有问题，回查 frame-0 距离后再上生产。
