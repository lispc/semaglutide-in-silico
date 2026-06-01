# exp-D: Linker 变体对 GLP-1R ECD 结合的影响

> 对应决策 4：γGlu-2×OEG 如何隔离脂链、恢复 GLP-1R 亲和力？
>
> 参考文献：Lau et al. (2015) J. Med. Chem. Table 3

## 科学问题

司美格鲁肽的 γGlu-2×OEG linker 将 C18 二酸脂链从肽骨架延伸出去，防止脂链干扰 GLP-1R ECD 结合界面。linker 太短（脂链贴靠 ECD）或太长（entropic penalty）都会降低 potency。

**核心假说**：linker 长度与 GLP-1R ECD 结合稳定性之间存在最优值（γGlu-2×OEG），与 Lau 2015 Table 3 的 EC50 趋势一致。

## 实验设计

5 个 linker 变体（均连接 C18 二酸，接在 Lys26）：

| # | Linker | Lau 2015 | hGLP-1R EC50 (pM) |
|---|--------|---------|-------------------|
| 1 | None | Cmpd 19 | 269 |
| 2 | γGlu | Cmpd 20 | 9.9 |
| 3 | γGlu-1×OEG | Cmpd 21 | 4.8 |
| 4 | γGlu-2×OEG | Semaglutide | 6.2 |
| 5 | γGlu-3×OEG | Cmpd 24 | 27.7 |

系统：GLP-1R ECD + Aib8,Arg34-GLP-1(10-35) + linker-C18 + TIP3P 水

## MD 协议

- 力场：ff14SB + GAFF2 + TIP3P, OpenMM 8.5.1
- 温度：310 K（Langevin, 1 ps⁻¹）
- 步长：2 fs, H 键约束（HBonds）
- 协议：minimization → NVT 0→100K (50ps) → NPT 100→310K (100ps) → NPT eq (200ps) → NPT production (500ns)
- 每条件 3 replica × 500 ns = 7.5 μs 总量
- 分析指标：ECD CA RMSD, C18-ECD 最短距离, linker RMSF

## 文件结构

```
exp-D/
├── README.md           # 本文件
├── exp-log.md          # 实验日志
├── tasks.md            # 任务跟踪
├── build/              # 分子构建脚本 + mol2 输入
│   ├── build_lya_variants.py  # 5 linker 变体批量构建
│   ├── strip_lya.py           # 剥离 LYA caps
│   ├── position_lnk.py        # Kabsch 定位到 Lys26 NZ
│   └── *_bcc.mol2             # AM1-BCC 电荷输出
├── tleap/              # tleap 拓扑构建脚本 + prmtop
│   ├── ecd_pep.pdb            # 清理后的 ECD+肽
│   ├── lya_link.frcmod        # 混合力场参数
│   └── lnk_*_pos.mol2         # 定位后的 linker mol2
├── md/                 # MD 生产
│   ├── run_md.py              # MD 脚本
│   └── {variant}/rep{N}/      # 轨迹 + 日志
└── analysis/           # 分析（待补充）
```

## 已知问题

1. **gglu NaN**：γGlu-only linker 的 C18 tail 碰撞 ECD 表面（PE=7e16），需调整初始旋转
2. **LNK 无共价键**：tleap bond 命令静默失败，LNK 原子靠非键相互作用维持（不影响 MD 稳定性）
3. **Arg34 突变**：Lys34→Arg 通过 tleap 模板补全，胍基几何为近似值

## 成功判据

- 无 linker 体系：ECD RMSD 最高，C18 tail 紧贴 ECD 表面
- γGlu-2×OEG：ECD RMSD 最低，C18 tail 远离 ECD
- 3×OEG：ECD 界面轻微不稳定（entropic penalty）
- 整体趋势与 Lau 2015 Table 3 EC50 数据定性一致

---

*最后更新：2026-06-01*
