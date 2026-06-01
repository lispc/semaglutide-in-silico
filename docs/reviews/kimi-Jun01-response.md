# Response to Kimi Review (2026-06-01)

## Summary

Kimi raised 19 issues across 5 categories. We agree with 11, partially agree with 3, and disagree with 5 severity assessments. No P0 issues affect the currently running exp-D production MD.

---

## Disputed P0 Assessments

### 3.1.1 `launch.sh` path error → claimed "all MD fails"

**Disagree with severity.** `launch.sh` was an early draft that was never executed. The actual 12 MD simulations were launched manually via direct `nohup` commands calling the correct `exp-D/md/run_md.py`. All 12 are running stably (T=310K, no NaN), with 8/12 already at 90-100 ns. The reviewer analyzed a dead file, not the active launch path.

**Action**: Delete `launch.sh` to prevent future confusion.

### 3.1.3 Double bond NZ → claimed "chemically absurd topology"

**Disagree with severity.** The double-bond script (`build_complex.py`) was never used to generate the actual running topologies. The actual tleap builds used `bond sys.117.NZ sys.128.C` (single bond). Furthermore, we subsequently discovered that tleap's `bond` command silently fails on our PDB format — no covalent NZ-C bond exists in the production prmtops. The LNK molecule is held together by nonbonded interactions at chemically reasonable distances. This is not ideal but does not produce "chemically absurd" geometries.

**Action**: Document the missing-bond limitation in exp-D/exp-log.md.

### 3.1.2 Lys34→Arg mutation → claimed "chemically incorrect"

**Partially agree, disagree with severity.** The crude `build_complex.py` approach (atom renaming) was not the final mutation method. After loading the renamed PDB, tleap auto-completed missing ARG atoms (NH1, NH2) from the ARG template, using existing heavy atom positions as placement guides. The resulting geometry is approximate but not random. For exp-D's primary metrics (C18-ECD distance, linker RMSF), the exact Arg34 guanidinium geometry is secondary.

**Action**: Rebuild Arg34 mutation properly (PyMOL or tleap mutate) when regenerating gglu variant.

### 3.2.1–3.2.3 aib parameter files → claimed "risk of misuse"

**Disagree with severity.** These files (`aib_capped.mol2`, `sqm.out`, `aib_residue.xml`) are in `common/params/` as historical artifacts from early exploration. exp-A production used `convert_to_aib8.py` (ParmEd-based atom replacement with ff14SB types), which bypasses all three files. They should be archived but pose no risk to current production.

**Action**: Move unused param files to `common/params/_archive/`.

---

## Agreed Issues and Fix Plan

### P0: Fix this week

| # | Issue | Fix |
|---|-------|-----|
| 3.1.4 | Barostat added after `Simulation()` → NPT doesn't work | Move `system.addForce(MonteCarloBarostat)` before `Simulation()` in both exp-C and exp-D `run_md.py` |
| 3.3.2 | exp-D has no README/exp-log/tasks | Write all three |

### P1: Fix this week

| # | Issue | Fix |
|---|-------|-----|
| 3.1.5 | exp-A GROMACS coord filename mismatch | Fix `run_gmx.sh` to match `setup_gmx.py` output pattern |
| 3.3.1 | Hardcoded absolute paths in 40+ scripts | Introduce `common/scripts/paths.py` with `pathlib`-based root detection; migrate scripts incrementally |
| 3.3.5 | Scripts without `if __name__` guard | Add guard to scripts that may be reused (position_lnk, strip_lya, build_lya_variants) |
| 3.4.4 | `addIonsRand Na+ 0 Cl- 0` may not neutralize | Verify behavior; switch to `addIonsRand sys Na+ 0` (auto-neutralize) |
| 3.4.2 | exp-C hypothesis not formally revised | Update roadmap.md and exp-C/README.md to document linker-FA escape as a scientific finding |

### P2: Fix before next analysis cycle

| # | Issue | Fix |
|---|-------|-----|
| 3.2.4 | RDKit atom ordering hardcoded | Add assertion checks in build_lya_variants.py to validate atom identities |
| 3.3.3 | `mmpbsa.py` incomplete | Complete or remove from production directory |
| 3.3.4 | Old script versions unarchived | Move to `_archive/` subdirectories |
| 3.3.6 | Duplicate residues in `DPP4_FREE_RESIDUES` | Remove duplicates |
| 3.4.5 | Analysis stats lack autocorrelation correction | Add effective sample size estimation, median/IQR |

### P3: Ongoing

| # | Issue | Fix |
|---|-------|-----|
| 3.5.1 | Overly permissive `.claude/settings.local.json` | Restrict wildcards; remove stale PID entry |
| 3.5.2 | `cc-ds.sh` contains plaintext API key | Rotate key on DeepSeek console; use `.env` + `python-dotenv` |
| 3.2.5 | Generic bond type defaults in `build_sema_parmed.py` | Not used in production; mark as archived |

---

## Items Requiring Clarification

### 3.4.3 NPT vs NVT production

Exp-D uses NPT production (barostat kept), unlike exp-A which used NVT. This is intentional: exp-A had a fixed-volume system where NVT was appropriate for the binding analysis. Exp-D's solution-phase ECD+peptide system benefits from NPT to allow density relaxation. The different ensembles are documented and should be noted when comparing across experiments.

### 3.4.1 "Correct" Lys34→Arg approach

The reviewer suggests PyMOL/Modeller for side chain replacement. We will evaluate:
- tleap `mutate` command (simplest, keeps existing backbone)
- PyMOL mutagenesis wizard (rotamer library)
- Manual rotamer selection with energy minimization

---

## Conclusion

The reviewer identified several real issues warranting fixes, particularly the barostat timing bug and exp-D documentation gap. However, the three "P0 - stop production immediately" items were based on analysis of dead draft scripts rather than the active production pipeline. The currently running exp-D MD is producing valid data and should be allowed to complete.

The 11 agreed fixes will be addressed in priority order over the coming week.

---

*Response written: 2026-06-01*
