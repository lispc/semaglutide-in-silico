# exp-F: GLP-1R 跨膜激活机制（全项目顶峰实验）

> 对应决策：决策 3+4（在完整膜受体环境中验证脂链和 Linker 的行为）
>
> **状态：Phase 2 — 膜系统构建 ✅ 完成（2026-06-06）**
>
> **前置依赖**：[x] exp-D 结论补齐或明确能力边界
> > - [x] 完整司美格鲁肽体系（受体 + 肽 + γGlu-2×OEG-C18 二酸）的构建管线验证
> > - [x] exp-C 的"linker-FA 逃逸"发现被整合进设计逻辑
>
> **战略定位**：
> exp-F 是全项目最重要的实验——它把故事放在正确的生物学上下文中（完整 GPCR + 膜 + Gs）。
>
> 但 review `docs/reviews/claude-Jun06.md` 指出一个关键问题：
> > "项目目前没有任何一个体系跑的是完整司美格鲁肽（exp-C 无受体、exp-D 无 HSA）。四步决策链里'脂链 vs 受体的空间竞争'这一核心矛盾，只有在完整分子（或 exp-F 的膜体系）里才能真正观察到。"
>
> 这意味着 exp-F 不只是"更贵的 exp-D"，而是**唯一能回答"司美格鲁肽为什么选 2×OEG"这个问题的实验**——因为 2×OEG 的真正优势体现在 BR ratio（HSA 结合 vs 受体结合的平衡），而这个平衡只有在同时包含受体和脂链的体系中才能观察到。

---

## Phase 历史

### Phase 0：构建管线验证（2026-06-06）
- 成功构建了完整司美格鲁肽 + GLP-1R 受体的溶剂化体系（363,332 atoms）
- 力场：ff14SB + GAFF2 + TIP3P
- 使用了 exp-D 验证的 LNK 参数化策略（loadMol2 + combine + remove N + bond NZ-C）
- 体系文件：`tleap/system.prmtop` / `system.inpcrd`

### Phase 1：稳定性验证（2026-06-06）
- Minimization + NVT→NPT equilibration (250 ps) 成功，无 NaN/能量爆炸
- 温度控制精准（100→310 K 各阶段 ±1 K）
- GPU 速度：~19 ns/day（RTX 3090, CUDA mixed）
- 最终势能：-4,145,861 kJ/mol（稳定）
- 文件：`md/rep1/equilibrated.pdb` / `equil.chk`

### Phase 2：膜嵌入系统构建（2026-06-06）
- **诊断 packmol-memgen 关键 bug**：`renumber=True` 导致脂质残基号与蛋白冲突，触发 `charmmlipid2amber.py` 的残基合并错误
- **手动 packmol 生成膜组件**：POPC:CHL1 (8:2) + 水 + 离子（562 lipids）
- **自定义 merge pipeline**：Python 脚本合并蛋白与膜，KD-tree 去除重叠脂质（移除 113，保留 449），分配唯一 chain ID，重编号残基
- **修复 PDB 列对齐 bug**：4-char resname 导致 chainID 移位，坐标列左偏 1 列
- **修复内部链断裂**：插入 5 处 TER 记录防止 tleap 生成超长 C-N 键
- **LNK 电荷归零**：将 +3.15842 归一化至 0.000000
- **最终体系**：312,501 atoms，449 lipids，79,615 waters，425 ions
- **tleap 构建成功**：`Errors = 0`，cpptraj 验证通过
- 文件：`membrane_build/system_final.{prmtop,inpcrd,pdb}`

---

## 建议的启动路径
1. ~~先构建一个"最小完整司美格鲁肽"体系（无膜，仅受体 + 肽 + linker + C18 二酸），验证构建管线~~ ✅ 已完成
2. ~~嵌入 POPC/胆固醇膜，修复 packmol-memgen 冲突，构建完整膜体系~~ ✅ 已完成
3. 跑 production MD：三组对比 WT (γGlu-2×OEG) / No-linker / Short-linker (仅 γGlu)

## 资源
- ~312k atoms/system（膜嵌入），OpenMM ~15-19 ns/day（RTX 3090）
- 100 ns × 3 replicas ≈ 16-21 天日历（1 GPU）

---

*创建日期: 2026-05-26*
*最后更新: 2026-06-06 — Phase 2 膜系统构建完成*
