# exp-D Final Results — Linker Variant Comparison with Significance Testing

Date: 2026-07-19. Script: `analysis/compare_linkers_v2.py` (supersedes `compare_linkers.py`, which had no significance testing and pointed at lost data).
Data: rebuilt topologies (`tleap/{variant}.prmtop`, 2026-07-17) + new trajectories (`md/{variant}/rep{1,2,3}/{variant}_traj.dcd`), 3 replicas × 100 ns per variant.

## Method

- **Metrics** (same definitions as v1 for cross-comparability):
  - CA RMSD — ECD CA only (resSeq 0–100), after alignment to ECD CA frame 0.
  - Tail-Prot — min distance from LNK tail carbons (lnk_c[2:], i.e. C14 onward; amide C11 and C13 excluded) to any protein CA.
- **dt = 0.1 ns** from DCDReporter interval (50,000 steps × 2 fs); DCD timestamps not read (unreliable).
- **Equilibration**: last 50% of each replica (frames 501–1002, 50–100 ns).
- **Frame-0 topology check**: NZ(Lys26)–C11(LNK amide) distance. All 12 replicas: 1.34–1.43 Å ✓ (fail threshold 3 Å; none skipped).
- **Statistics** (`common/lib/stats.py`, first real use of `correlated_t_test`):
  - Primary test: correlated t-test on concatenated per-frame equilibrated series (3 replicas concatenated; replicas are independent runs, pairing rep *i* of A with rep *i* of B is arbitrary — the test is used as a mean-difference test with autocorrelation-corrected n_eff).
  - Control: Welch t-test on per-replica means (n=3 vs 3 — low power, reported for reference only).
  - Multiplicity: 3 primary comparisons; Bonferroni α = 0.0167 noted alongside raw α = 0.05.

## Per-variant summary (equilibrated last 50 ns)

mean±SD over replica means; median/IQR pooled over equilibrated frames; n_eff = Σ per-replica n_eff (autocorrelation-corrected).

| Variant | Lau EC50 (pM) | Tail-Prot (Å) | median | IQR | n_eff | rep-CV | CA RMSD (Å) | median | IQR | n_eff | rep-CV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| No linker (19) | 269 | **5.04±0.26** | 4.50 | 1.90 | 48 | 5.1% | 1.77±0.02 | 1.74 | 0.28 | 57 | 1.4% |
| γGlu-1×OEG (21) | 4.8 | 3.74±0.14 | 3.68 | 0.54 | 45 | 3.6% | 1.71±0.04 | 1.70 | 0.20 | 85 | 2.2% |
| γGlu-2×OEG (Sema) | 6.2 | 3.95±0.07 | 3.93 | 0.30 | 508 | 1.8% | 1.68±0.24 | 1.65 | 0.39 | 101 | 14.1% |
| γGlu-3×OEG (24) | 27.7 | 3.65±0.12 | 3.66 | 0.30 | 97 | 3.4% | 1.82±0.14 | 1.81 | 0.27 | 94 | 7.5% |

## Primary question: gglu_2oeg (semaglutide linker) vs each other variant — Tail-Prot

Correlated t-test on concatenated per-frame series:

| Comparison | Δ (Å) | t | df | p | n_eff | Welch (rep means, n=3) t (df) p |
|---|---:|---:|---:|---:|---:|---|
| 2×OEG vs no_linker | −1.09 | −4.95 | 41.6 | **1.3×10⁻⁵** | 43 | −7.07 (2.3) p=0.013 |
| 2×OEG vs gglu_1oeg | +0.21 | +2.58 | 32.1 | **0.0147** | 33 | +2.39 (3.0) p=0.097 |
| 2×OEG vs gglu_3oeg | +0.30 | +3.34 | 15.2 | **0.0044** | 16 | +3.66 (3.2) p=0.032 |

All three comparisons are significant at α=0.05, and vs no_linker / vs 3×OEG also pass Bonferroni α=0.0167 (vs 1×OEG: p=0.0147 < 0.0167, marginal pass).

## Secondary: gglu_3oeg vs others — CA RMSD

| Comparison | Δ (Å) | t | df | p | n_eff | Welch (rep means) t (df) p |
|---|---:|---:|---:|---:|---:|---|
| 3×OEG vs no_linker | +0.05 | 0.85 | 20.2 | 0.41 | 21 | 0.67 (2.1) p=0.57 |
| 3×OEG vs gglu_1oeg | +0.11 | 1.54 | 13.6 | 0.15 | 15 | 1.35 (2.3) p=0.29 |
| 3×OEG vs gglu_2oeg | +0.14 | 1.00 | 5.1 | 0.36 | 6 | 0.92 (3.2) p=0.42 |

No significant CA-RMSD difference anywhere (all p > 0.14). Note gglu_2oeg CA RMSD replica CV = 14% (rep1 1.95 vs rep3 1.50 Å) — ECD flexibility varies between replicas at this timescale.

## Conclusions

1. **2×OEG differs significantly from every other variant in Tail-Prot distance** — this is the first time this claim is backed by a test. But the *pattern* is not "2×OEG is special":
   - vs no_linker: **−1.09 Å, p≈1×10⁻⁵** — the dominant, clearly significant effect. Any γGlu-(OEG)n linker holds the C18 tail ~1 Å closer to the ECD surface than the bare C18.
   - vs 1×OEG / 3×OEG: +0.21 / +0.30 Å, p = 0.015 / 0.004 — statistically detectable thanks to n_eff-corrected testing, but the effect is **0.2–0.3 Å, i.e. sub-Å and close to thermal fluctuation scale**; physical relevance is doubtful, and the sign says 2×OEG sits slightly *farther* from the protein than 1×/3×OEG.
2. **The historical preliminary result does not replicate.** Old (lost) trajectories suggested 2×OEG uniquely compact (3.8 Å vs 4.1–4.6 Å, within error bars, untested). New data: ordering is 3×OEG (3.65) < 1×OEG (3.74) < 2×OEG (3.95) ≪ no_linker (5.04). 2×OEG is intermediate, not the most compact.
3. **Tail-Prot distance does not track Lau 2015 Table 3 potency** (1×OEG ≈ 2×OEG > 3×OEG ≫ no linker). The only qualitatively consistent extreme is no_linker: farthest tail *and* 40–50× least potent (269 pM). But 3×OEG is the *closest* in simulation while being 4.5× *less* potent than 2×OEG — so "closer tail = more potent" is falsified as a mechanism for discriminating 1×/2×/3×OEG. If linker length matters for potency, it is not via mean tail–ECD proximity at the 100-ns timescale (alternative: tail–membrane/HSA partitioning, or rare conformational events not sampled here).
4. **The old "3×OEG highest CA RMSD → entropic penalty" claim is not supported** either: 3×OEG ECD RMSD (1.82 Å) is statistically indistinguishable from all other variants (all p > 0.14).

**One-sentence answer:** yes, γGlu-2×OEG is statistically distinguishable from the other variants in Tail-Prot distance (p ≤ 0.015 for all three pairwise tests), but the only large and clearly meaningful difference is "any linker vs no linker" (−1.1 Å); among the three linker lengths the differences are 0.2–0.3 Å, opposite in sign to the old preliminary claim, and do not order the way Lau 2015 potency does.

## Reproducibility

- Raw per-frame series: `analysis/per_frame_series.npz` (`{variant}_r{rep}_{ca_rmsd,tail_prot}`).
- All numbers: `analysis/test_results.json`; console log: `analysis/compare_linkers_v2.log`.
- Figure: `analysis/tail_prot_compare.png` (Tail-Prot distribution, variant × replica facets).

## Caveats

- n=3 replicas per variant; per-replica-mean Welch tests have df≈2–3 and are underpowered — they corroborate vs no_linker but miss the small 1×/3×OEG effects.
- Correlated t-test pairing of concatenated independent replicas is arbitrary; conclusions rest on the n_eff-corrected mean difference, not on genuine pairing.
- LNK partial charges are per-element placeholder averages (carried over from the original build, identical across variants), not AM1-BCC — acceptable for comparative ranking, not for absolute energetics.
- 100 ns × 3 may undersample slow tail–surface unbinding (no_linker rep1/rep2 show excursions to 6–9 Å, see figure).
