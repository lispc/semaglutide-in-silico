# exp-A 实验日志

> 只追加，不删除。每次记录含日期、时间、操作内容和结果。

---

## 2026-05-26 — 实验目录初始化

- 创建 `exps/exp-A/` 目录结构
- 编写 README.md（实验设计）、tasks.md（任务清单）、exp-log.md（本文件）
- **状态**：设计阶段，等待用户确认后启动 Phase 0

## 2026-05-26 — 设计调整（用户反馈）

调整项：
- **肽长度**：9 残基 (7-15) → 全长 31 残基 (7-37)。N 端 backbone restraint (res 1-6, 100 kJ/mol/nm²) 防长肽漂移
- **温度**：300 K → **310 K**（生理温度）
- **力场**：ff19SB → **ff14SB**（ff19SB 的残基特异性 CMAP 类型不支持非标准 Aib 残基）
- **引擎**：纯 OpenMM，不用 GROMACS

## 2026-05-26 — tleap 拓扑构建

- WT: tleap → prmtop/inpcrd → OpenMM 最小化。PE: 1.81×10¹⁵ → -1.73×10⁶ kJ/mol。✓
- Aib8: tleap 无法正确处理 AIB mol2（mol2 无 IC 表，H 原子坐标为 0）。
  改用 ParmEd 从 WT 拓扑修改 ALA→AIB：去 HA、加 CB2+3H、调电荷。无 clash。OpenMM 最小化成功。
- 服务器：tleap (AmberTools 24.8) from env `cgas-md`，OpenMM 8.5.1 from env `gmx`

## 2026-05-26 — 200 ns 生产 MD 启动

- **15:30** WT (GPU 0) 和 Aib8 (GPU 1) OpenMM 同时启动
- **17:00** WT (GPU 2) 和 Aib8 (GPU 3) GROMACS 交叉验证启动
- OpenMM 脚本: `common/scripts/run_md.py`
- GROMACS 脚本: `exps/exp-A/gmx/run_gmx.sh`
- OpenMM 约束: DPP-4 CA 弱约束 + 肽 N 端 BB 约束
- GROMACS: 无约束（仅用作引擎交叉验证）
- 310 K, NVT production, 2 fs step, checkpoint 每 1 ns

---

*维护者：Claude Code*
*最后更新：2026-05-26*
