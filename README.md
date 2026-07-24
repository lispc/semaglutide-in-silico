# Semaglutide in Silico: Reconstructing the Rational Design of a Once-Weekly GLP-1 Analogue

[![GitHub](https://img.shields.io/badge/repo-semaglutide--in--silico-blue)](https://github.com/lispc/semaglutide-in-silico)

An atomistic, computational re-creation of the four-step design logic that transformed native GLP-1 into semaglutide, a once-weekly glucagon-like peptide-1 (GLP-1) analogue used for type 2 diabetes and obesity.

---

## Table of Contents

- [Project Goal](#project-goal)
- [The Four Design Decisions](#the-four-design-decisions)
- [Repository Structure](#repository-structure)
- [Experiments](#experiments)
- [Current Status](#current-status)
- [Technical Stack](#technical-stack)
- [Key Findings So Far](#key-findings-so-far)
- [Known Limitations](#known-limitations)
- [Reproduction & Usage](#reproduction--usage)
- [References](#references)
- [License](#license)

---

## Project Goal

Semaglutide (Ozempic, Wegovy, Rybelsus) is the result of a tightly coupled chain of rational mutations and chemical modifications. The goal of this project is to reproduce that chain **at atomic resolution** using modern molecular dynamics (MD) simulations, free-energy methods, and cross-validation practices.

We are not only trying to explain *what* the final molecule looks like, but *why each preceding modification was necessary* — and how the four design decisions interact.

Core references:

- Lau, J. *et al.* (2015) *J. Med. Chem.* **58**, 7370–7380 — original SAR paper.
- Knudsen, L. B. & Lau, J. (2019) *Front. Endocrinol.* **10**, 155 — design narrative review.

---

## The Four Design Decisions

The project follows the causal chain used by Novo Nordisk to engineer semaglutide from native GLP-1:

```
native GLP-1
    │ t½ ~1–2 min; cleaved by DPP-4 at Ala8-Glu9
    ▼ Decision 1: Ala8 → Aib8 (steric block of DPP-4)
Aib8-GLP-1
    │ still cleared rapidly by kidneys
    ▼ Decision 2: Lys34 → Arg34 (site-specific acylation at Lys26)
Aib8,Arg34-GLP-1(Lys26-acyl)
    │ albumin binding extends t½, but…
    │   • C16 mono-acid affinity is too low
    │   • long mono-acids aggregate and are poorly soluble
    │   • bulky acyl chain clashes with the GLP-1 receptor
    ▼ Decision 3: C16 mono-acid → C18 di-acid
    ▼ Decision 4: γGlu-2×OEG flexible linker (decouples acyl chain from receptor)
Semaglutide
    │ t½ ~165 h; once-weekly dosing
```

Each decision is tested in a dedicated experiment (A–F) with MD-based structural evidence.

---

## Repository Structure

```
semaglutide-in-silico/
├── common/
│   ├── lib/              # Shared Python utilities (stats.py, etc.)
│   ├── params/           # Force-field parameters (Aib, linker, acyl chain)
│   └── scripts/          # Reusable MD/build wrappers
├── docs/
│   ├── best-practice.md  # Lessons learned from earlier projects
│   ├── reviews/          # External / self code & doc reviews
│   ├── report_draft.md   # Manuscript outline
│   └── lessons_learned.md
├── exps/
│   ├── exp-A/            # DPP-4 binding / Aib8 steric block
│   ├── exp-B/            # Lys34 → Arg34 acylation control
│   ├── exp-C/            # HSA acyl-chain binding optimization
│   ├── exp-D/            # Linker length / OEG repeat effects
│   ├── exp-E/            # Full SAR-table in-silico reproduction
│   └── exp-F/            # Full-length GLP-1R in membrane
├── roadmap.md            # Detailed Chinese-language project roadmap
└── README.md             # This file
```

---

## Experiments

| Experiment | Decision | Question | System | Method |
|------------|----------|----------|--------|--------|
| **A** | 1 | Why does Aib8 block DPP-4? | DPP-4 + GLP-1(7-37) WT vs Aib8 | cMD + MM-PBSA |
| **B** | 2 | Why must it be K34R? What if Lys34 were acylated? | GLP-1R ECD + WT / R34 / acyl-K34 | cMD |
| **C** | 3 | Why is C18 di-acid optimal for HSA binding? | HSA FA3 site + acyl-chain variants | cMD + MM-GBSA |
| **D** | 4 | How does γGlu-2×OEG resolve steric clashes? | GLP-1R ECD + linker variants | cMD |
| **E** | All | Can the pipeline quantitatively reproduce Lau Tables 1–4? | ~15–20 analogues | cMD + correlation analysis |
| **F** | 3 + 4 | How do acyl chain and linker behave in the full membrane receptor? | GLP-1R-TMD + Gs + POPC/CHOL | cMD |

Planned aggregate sampling: ~47 µs across 35+ systems (4× RTX 3090).

---

## Current Status

Last major update: 2026-07-24

| Item | Status | Notes |
|------|--------|-------|
| Phase 0 infrastructure | ✅ Complete | ff14SB + GAFF2 + TIP3P pipeline; Aib parameterization; GROMACS/OpenMM cross-checks |
| exp-A (DPP-4) | ✅ Complete (1 replica) | Aib8 double-methyl consistently pushes the peptide away from the DPP-4 active site; formal criterion (WT ≤3.5 Å) not met — conclusions are relative. MM-PBSA done (2026-07-19): global ΔΔG opposite to expected sign (artifact-dominated, recorded as a negative result); P2-site decomposition supports the local desolvation penalty (+3.1/+3.4 kcal/mol) |
| exp-C (HSA acyl chain) | ✅ Complete + MM-GBSA + chain series | Free FAs anchor FA3 via distal-carboxyl double salt bridge (ARG346+ARG483, matching Liu 2025); MM-GBSA (2026-07-20): di-acid beats mono-acid (ΔΔG = −10.8 ± 2.6 GB / −18.5 ± 4.8 PB kcal/mol, p≈0.03) via second salt bridge to ARG408; attached linker still drives tail escape (hydrophilic OEG). Chain-length series (2026-07-24): C12–C22 di-acids ×3 replicas + MM-GBSA U-curve — C18 sits in the top cluster (GB −87.4) but is not uniquely optimal (C22 strongest at −93.2, within method error); C16 di-acid is a consistent local minimum; FAH residue contribution grows monotonically with chain length |
| exp-G (HSA/receptor competition) | 🔄 Running (2/3 replicas) | First ternary system (HSA + semaglutide + ECD, ~142k atoms, FA3 pre-anchored): pilot 20 ns shows the linker accelerates FA3 undocking ~100× (~2 ns), yet the tail stays on the HSA surface (2.53 ± 0.14 Å) while the peptide stays locked on ECD (1.34 ± 0.03 Å) — direct evidence that γGlu-2×OEG lets one molecule hold albumin and receptor simultaneously |
| exp-D (linker) | ✅ Complete + stats + MM-GBSA | First statistical test (correlated t-test, 2026-07-19): linker-vs-no-linker effect significant (−1.1 Å, p≈1e-5); OEG-length differences only 0.2–0.3 Å and do not order like Lau Table 3 potency — old preliminary claims not replicated. MM-GBSA (2026-07-20) confirms: only the linker/no-linker contrast is robust; length ordering mismatches potency; Spearman with EC50 positive but underpowered (n=4) |
| exp-F membrane system | ✅ Complete (1 replica) | GROMACS 200 ns + OpenMM 224 ns cross-validation; membrane stabilizes the receptor (protein CA RMSD 8.7 / 4.7 Å per engine vs ~32 Å without membrane); acyl tail stays in water (40.6–42.6 Å from membrane center, zero lipid contacts) |
| exp-F minimal ECD model | ✅ Complete (1 replica) | 99.3 ns; internal RMSDs stable (ECD 1.4 Å, peptide 2.7 Å), but HSA never approaches the acyl tail (COM ~109 Å) — HSA/receptor competition not yet observed |
| exp-B (K34R) | ⏸ Paused | ECD-peptide complexes separate in all crystal structures; waiting for AF3/Boltz-1 alternative |
| exp-E (full SAR) | ⏳ Not started | Pending Phase 1/2 validation |

See `roadmap.md` (Chinese) and `docs/reviews/` for the most detailed status and methodological boundaries.

---

## Technical Stack

| Component | Choice | Rationale |
|-----------|--------|-----------|
| Primary MD engine | **OpenMM 8.5.1** (small/medium systems); **GROMACS 2026** (large membrane system) | Engine speed is system-dependent: OpenMM ~178 ns/day on the ~140k-atom exp-A system, but GROMACS is 2.3× faster on the 312k-atom membrane system (96 vs 42 ns/day — see best-practice-v2 #27) |
| Cross-validation engine | **GROMACS** / **OpenMM** (mutual) | Mature analysis toolchain; dual-engine sanity checks |
| Protein force field | **ff14SB** | ff19SB CMAP types are incompatible with non-standard Aib |
| Water model | **TIP3P** | Consistent with ff14SB |
| Non-standard residues | **AmberTools tleap + ParmEd + GAFF2** | Manual mol2 construction because antechamber (`sqm`/`bondtype`) fails for linker/acyl groups |
| Structure prediction | **AF3, Boltz-2, Chai-1** | Cross-validation of low-confidence interfaces |
| Docking | **LightDock, Rosetta, ClusPro** | Multi-method pose consensus |
| Analysis | **MDAnalysis, mdtraj, cpptraj, PyMOL** | Independent tool cross-checks |
| Statistics | `common/lib/stats.py` | Autocorrelation-corrected summaries and correlated t-tests |

Production protocols:

- Temperature: 310 K
- Pressure: 1 bar (Monte Carlo barostat)
- Timestep: 2 fs with SHAKE on H-bonds
- Cutoff: 10 Å with PME for long-range electrostatics
- Replica count: ≥3 per condition

---

## Key Findings So Far

1. **Aib8 steric block (exp-A)** — The Aib8 double-methyl pushes the peptide away from the DPP-4 active site on every geometric metric (original run: catalytic distance +1.0 Å, S1-pocket contacts −49, CB→Trp629 +1.7 Å; 2026-07 re-run reproduces the direction: catalytic distance 3.70 vs 5.25 Å WT/Aib8). MM-PBSA on the re-run (2026-07-19) does **not** show globally weaker binding — the global ΔΔG is artifact-dominated and recorded as a negative result: DPP-4 resistance is interpreted as catalytic-pose exclusion, not global affinity loss. The P2-site decomposition is unfavorable for Aib8 (+3.1/+3.4 kcal/mol), as hypothesized.

2. **HSA acyl-chain behavior (exp-C)** — Free fatty acids anchor HSA FA3 via a distal-carboxyl double salt bridge (ARG346+ARG483, 2.72–2.75 Å, matching Liu 2025's R348/R485 picture); MM-GBSA (2026-07-20) shows the C18 di-acid binds significantly stronger than the mono-acid (ΔΔG = −10.8 ± 2.6 GB / −18.5 ± 4.8 PB kcal/mol, p≈0.03, consistent across 3/3 replicas and both methods), driven by a second salt bridge from the proximal carboxyl to ARG408 — directly supporting the mono→di-acid design decision. Once the linker is attached, the fatty-acid tail escapes 32–41 Å from the protein (hydrophilic OEG is the dominant driver), implying that linker design controls both albumin anchoring and receptor accessibility.

3. **Membrane environment is essential (exp-F)** — In a solvent-only system the GLP-1R transmembrane domain drifts dramatically (RMSD ~32 Å). Embedding the receptor in a POPC/cholesterol bilayer reduces protein CA RMSD to 8.7 Å (GROMACS 200 ns) / 4.7 Å (OpenMM 224 ns) and keeps the peptide N-terminus near the receptor.

4. **C18 di-acid does not insert into the membrane** — In the membrane-bound receptor system the acyl tail stays in the aqueous phase, 40.6–42.6 Å from the membrane center with zero lipid contacts (consistent across both engines), matching the exp-C "linker-FA escape" observation.

5. **Statistical rigor is being added incrementally** — `common/lib/stats.py` provides Geyer IACT-based effective sample sizes, replica CV and correlated t-tests; n_eff/CV reporting is wired into the exp-A/exp-D analyses, but formal t-tests between conditions have not yet been run.

6. **Chain-length optimality is a plateau, not a sharp optimum (exp-C chain series, 2026-07-24)** — MM-GBSA over C12–C22 di-acids (3 replicas each, true carbon counts after the naming erratum) shows binding strengthens from C12 (−71.3 GB kcal/mol) into a top cluster spanning C18–C22 (−87.4 … −93.2), with C16 a consistent local minimum in both GB and PB; per-residue decomposition attributes the trend to monotonically growing FAH-pocket contributions (−13.4 → −29.6 kcal/mol). Semaglutide's C18 di-acid is therefore "safely inside the optimal cluster" rather than a unique energy minimum — the choice plausibly also reflects synthesis, solubility and receptor-activity constraints. Di-acid beats mono-acid at equal length (C16: −76.2 vs −70.5 GB), and the γGlu-2×OEG linker lets the ternary complex (exp-G pilot) keep the fatty-acid tail on the HSA surface while the peptide stays bound to the receptor ECD.

---

## Known Limitations

The project is intentionally transparent about what has **not** yet been demonstrated:

- **The full competition system exists only as exp-G (running).** exp-C lacks receptor; exp-D lacks HSA; exp-F lacks HSA. exp-G (HSA + peptide-linker-FA + ECD, FA3 pre-anchored) closes this gap — pilot + 2/3 production replicas done as of 2026-07-24 — but its residence-partition analysis is still pending, so the competition is not yet quantified across replicas.
- **Sampling is below the original roadmap target** for exp-C and exp-D (100 ns achieved vs 300–500 ns planned). Slow degrees of freedom may not be converged; conclusions are labelled preliminary.
- **MM-GBSA is now delivered for exp-A/C/D** (exp-A: negative result; exp-C: di-acid > mono-acid, ΔΔG = −10.8/−18.5 kcal/mol; exp-D: only the linker/no-linker contrast is robust, consistent with the t-test). **Chain-length naming erratum (2026-07-20)**: the legacy "C18" free-FA controls are actually C19 monoacid / C20 diacid (original scripts named chains by CH₂ count); the mono-vs-diacid conclusion is unaffected. The follow-up chain-length series (true C12/C14/C16/C18/C20/C22 di-acids + C16 monoacid, 3 replicas each) completed MD on 2026-07-24 with MM-GBSA U-curve delivered (`exps/exp-C/analysis/mmgbsa/chain-series/`); the final c18true replica 3 point is still being computed.
- **Ensemble inconsistency**: exp-A used NVT production; exp-C/D used NPT production. Cross-experiment energy comparisons must account for ~2% density differences.
- **Statistical pipeline is still being wired in**: early analyses reported mean±std without autocorrelation correction; this is being retrofitted via `common/lib/stats.py`.

See `roadmap.md` §“能力边界与方法学一致性” and `docs/reviews/claude-Jun06.md` for the full methodological audit.

---

## Reproduction & Usage

### Environment

Conda/Miniforge environments used in this project include:

- `gmx` — OpenMM 8.5.1, MDAnalysis, mdtraj, numpy, scipy
- `cgas-md` — AmberTools 24.8 (`tleap`, `antechamber`, `parmchk2`)
- `rosetta` — Rosetta 2026.15
- `boltz` — Boltz-2 structure prediction

### Quick-start: build a minimal system

```bash
# Example: exp-A Aib8 peptide + DPP-4 build flow
cd exps/exp-A/tleap
python convert_to_aib8.py          # mutate Ala8 to Aib8 in the GLP-1 sequence
python fix_aib8_prmtop.py          # patch Aib parameters into topology

cd ../md
python run_md.py --system aib8     # production MD (after minimization/equilibration)
```

### Analysis

```python
from common.lib.stats import summarize, correlated_t_test
import numpy as np

# Load a 1-D observable (e.g., catalytic distance)
d = np.loadtxt("distances.dat")
print(summarize(d, name="Ser630-OG_to_Aib8-C"))
```

### Notes for Reproducers

- Always perform **frame-0 validation**: clash check, key distances, and initial energy.
- GROMACS trajectories must be unwrapped with `gmx trjconv -pbc whole` before analysis.
- Check `docs/best-practice.md` for hard-won lessons on multi-engine validation, PDB column alignment, non-standard residue parameterization, and statistical testing.

---

## References

1. Lau, J. *et al.* Discovery of the Once-Weekly Glucagon-Like Peptide-1 (GLP-1) Analogue Semaglutide. *J. Med. Chem.* **2015**, 58, 7370–7380.
2. Knudsen, L. B. & Lau, J. The Discovery and Development of Liraglutide and Semaglutide. *Front. Endocrinol.* **2019**, 10, 155.
3. Frimann, D. A. *et al.* Molecular Dynamics Insights into the Binding Interactions of Acylated GLP-1 Analogues with the GLP-1 Receptor. *J. Biomol. Struct. Dyn.* **2023**, 41(11), 5007–5021.
4. Liu, Y. *et al.* Molecular Dynamics Study of Semaglutide Binding to Human Serum Albumin. *J. Biomol. Struct. Dyn.* **2025**.
5. Sønderby, P. *PhD Thesis* (2016) — HSA self-interaction and liraglutide-HSA binding.

---

## License

This is a research project. Code and data are provided for academic and educational use. Please cite the original Novo Nordisk papers (Lau 2015; Knudsen & Lau 2019) if you use the scientific rationale or structural hypotheses presented here.

---

*Project maintained by the semaglutide-in-silico team. For the detailed Chinese-language roadmap, see `roadmap.md`.*
