# exp-C 任务清单

> 状态: □ 未开始 | ● 进行中 | ✓ 完成 | ✗ 跳过

---

## Phase 0: 游离 FA 验证（已完成，保留为基线）

- [✓] **0.0** HSA 1E7G 下载 + FA3 位点识别 (MYR 1003)
- [✓] **0.1** C18 monoacid 参数化 + tleap 拓扑
- [✓] **0.2** C18 monoacid ×3 replica × 100 ns → 远端羧基锚定 ARG482 2.8 Å
- [✓] **0.3** C18 diacid 参数化修复 (O1P-O2P 1.2→2.2 Å) + tleap 拓扑
- [✓] **0.4** C18 diacid ×3 replica × 100 ns → 仅单点锚定，近端羧基游离
- [✓] **0.5** 文献对照：确认 linker 对双点锚定必要性

## Phase 1: Linker-FA 复合体参数化

- [ ] **1.0** 构建 γGlu-2×OEG linker 3D 结构
- [ ] **1.1** 构建 linker-C18 二酸复合体 (γαE2-C18 diacid)
- [ ] **1.2** 构建 linker-C16 二酸、linker-C20 二酸
- [ ] **1.3** 构建 linker-C16 单酸、linker-C18 单酸
- [ ] **1.4** GAFF2 参数化（mol2 + 电荷归一化）
- [ ] **1.5** 各变体定位至 FA3（以 MYR 1003 为模板）

## Phase 2: Linker-C18 验证 MD

- [ ] **2.0** tleap 构建 HSA + linker-C18 二酸拓扑
- [ ] **2.1** linker-C18 二酸 ×3 replica × 100 ns
- [ ] **2.2** 初步分析：双点锚定是否成立？vs 游离 FA 对比

## Phase 3: 全 Linker-FA 集合 MD

- [ ] **3.0** linker-C12 二酸 ×3 replica × 300 ns
- [ ] **3.1** linker-C14 二酸 ×3 replica × 300 ns
- [ ] **3.2** linker-C16 二酸 ×3 replica × 300 ns
- [ ] **3.3** linker-C18 二酸 ×3 replica × 300 ns（阳性对照）
- [ ] **3.4** linker-C20 二酸 ×3 replica × 300 ns
- [ ] **3.5** linker-C16 单酸 ×3 replica × 300 ns
- [ ] **3.6** linker-C18 单酸 ×3 replica × 300 ns

## Phase 4: 全量分析

- [ ] **4.0** MM-GBSA 结合自由能（所有体系，最后 100 ns）
- [ ] **4.1** ΔG_bind vs 链长曲线（核心图表）
- [ ] **4.2** 双羧基锚定 occupancy 定量
- [ ] **4.3** Per-residue 能量分解（linker vs FA vs 远端羧基）
- [ ] **4.4** 与 Liu 2025 的 -75.33 kcal/mol 方法学交叉验证
- [ ] **4.5** 与游离 FA 基线对比（linker 的定位贡献）
- [ ] **4.6** 综合图表

## 已知风险

| 风险 | 缓解 |
|------|------|
| Linker-FA 参数化复杂 (~50 atoms) | 使用标准 GAFF2 官能团电荷，参考 Frimann 2023 参数化方案 |
| Linker 柔性导致收敛慢 | 延长至 300 ns，取最后 100 ns 分析 |
| 近端羧基可能仍无法锚定 | 若 linker 长度不够到达第二正电荷残基，则验证 linker 长度效应 |

---

*创建日期：2026-05-27*
*最后更新：2026-05-28 — 升级为 linker-脂链复合体*
