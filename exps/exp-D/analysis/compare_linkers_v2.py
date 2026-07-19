#!/usr/bin/env python3
"""exp-D final analysis: 4 linker variants, 3 replicas x 100 ns, WITH significance testing.

Differences vs compare_linkers.py (2026-06-06):
- New rebuilt topologies (tleap/{variant}.prmtop) and new trajectories (md/{variant}/repN/).
- Frame-0 topology/trajectory match check: NZ(Lys26)-C11(LNK amide) must be ~1.5 A;
  replicas with >3 A are reported and skipped.
- CA RMSD on ECD CA only (resSeq 0-100), aligned to ECD CA frame 0.
- Tail-Prot = min distance from LNK tail carbons (lnk_c[2:], i.e. C14+; C11/C13 skipped)
  to any protein CA — identical definition to compare_linkers.py for cross-comparability.
- dt = 0.1 ns derived from DCDReporter interval (50000 steps x 2 fs = 100 ps);
  DCD timestamps are unreliable and are NOT read.
- Last 50% of each replica treated as equilibrated (same convention as v1).
- Significance: common.lib.stats.correlated_t_test on concatenated per-frame
  equilibrated series (replicas are independent runs; pairing rep i of A with rep i
  of B is arbitrary — noted in RESULTS.md). Control: unpaired t-test on per-replica
  means (n=3, low power).
"""
import mdtraj as md, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os, sys, json

REPO_ROOT = "/home/scroll/zzhang/semaglutide-in-silico"
sys.path.insert(0, f"{REPO_ROOT}/common/lib")
import stats  # summarize, replica_cv, correlated_t_test

EXP_D = f"{REPO_ROOT}/exps/exp-D"
TLEAP = f"{EXP_D}/tleap"
MD_DIR = f"{EXP_D}/md"
OUT = f"{EXP_D}/analysis"
os.makedirs(OUT, exist_ok=True)

VARIANTS = ["no_linker", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]
LABELS = {"no_linker": "No linker (19)", "gglu_1oeg": "γGlu-1×OEG (21)",
          "gglu_2oeg": "γGlu-2×OEG (Sema)", "gglu_3oeg": "γGlu-3×OEG (24)"}
COLORS = {"no_linker": "#E53935", "gglu_1oeg": "#FB8C00",
          "gglu_2oeg": "#43A047", "gglu_3oeg": "#1E88E5"}
LAU_EC50 = {"no_linker": 269, "gglu_1oeg": 4.8, "gglu_2oeg": 6.2, "gglu_3oeg": 27.7}  # pM, Lau 2015 Table 3

DT_NS = 0.1  # 50000 steps * 2 fs
ECD_MAX_RESSEQ = 100  # residues 0..100 = ECD; 101..125 = peptide; 126 = LNK
NZC_FAIL_A = 3.0      # frame-0 NZ-C11 above this => topology/trajectory mismatch

results = {v: {} for v in VARIANTS}
skipped = []

for vname in VARIANTS:
    prmtop = f"{TLEAP}/{vname}.prmtop"
    for rep in [1, 2, 3]:
        dcd = f"{MD_DIR}/{vname}/rep{rep}/{vname}_traj.dcd"
        if not os.path.exists(dcd):
            print(f"MISSING: {dcd}"); skipped.append((vname, rep, "missing dcd")); continue
        print(f"Loading {vname} rep{rep} ...", flush=True)
        t = md.load(dcd, top=prmtop)
        print(f"  {t.n_frames} frames (dt={DT_NS} ns from reporter interval)")

        # --- selections ---
        prot_ca = t.topology.select("protein and name CA")
        ecd_ca = np.array([a.index for a in t.topology.atoms
                           if a.residue.is_protein and a.name == 'CA' and a.residue.resSeq <= ECD_MAX_RESSEQ])
        lnk_c = t.topology.select("resname LNK and element C")
        lnk_tail = lnk_c[2:]  # skip amide C11 + C13 (same as v1)
        nz_atoms = [a.index for a in t.topology.atoms
                    if a.residue.name == 'LYS' and a.name == 'NZ' and a.residue.resSeq == 116]
        assert len(nz_atoms) == 1 and len(lnk_c) >= 3 and len(ecd_ca) > 0
        nz, lnk_c0 = nz_atoms[0], lnk_c[0]  # C11 amide carbon

        # --- frame-0 NZ-C11 check ---
        nz_c0 = np.linalg.norm(t.xyz[0, nz] - t.xyz[0, lnk_c0]) * 10
        print(f"  frame-0 NZ(Lys26)-C11 = {nz_c0:.2f} Å")
        if nz_c0 > NZC_FAIL_A:
            print(f"  FAIL: NZ-C11 > {NZC_FAIL_A} Å — topology/trajectory mismatch, skipping")
            skipped.append((vname, rep, f"NZ-C11={nz_c0:.2f} A")); continue

        # --- metrics (aligned on ECD CA) ---
        t.superpose(t[0], atom_indices=ecd_ca)
        ca_rmsd = md.rmsd(t, t[0], atom_indices=ecd_ca) * 10
        prot_xyz = t.xyz[:, prot_ca]
        tail_xyz = t.xyz[:, lnk_tail]
        tail_prot = np.sqrt(((tail_xyz[:, :, None, :] - prot_xyz[:, None, :, :]) ** 2).sum(-1)).min(axis=(1, 2)) * 10

        h = t.n_frames // 2  # last 50% = equilibrated
        results[vname][rep] = {'ca_rmsd': ca_rmsd, 'tail_prot': tail_prot,
                               'half': h, 'nframes': t.n_frames, 'nz_c0': nz_c0}
        for arr, lab in [(ca_rmsd[h:], 'CA_RMSD'), (tail_prot[h:], 'Tail-Prot')]:
            print(f"  {stats.format_summary(stats.summarize(arr, name=lab))}", flush=True)
        del t

if skipped:
    print("\nSKIPPED:", skipped)

# --- persist raw per-frame series for reproducibility ---
np.savez(f"{OUT}/per_frame_series.npz",
         **{f"{v}_r{r}_{m}": results[v][r][m] for v in VARIANTS for r in results[v]
            for m in ('ca_rmsd', 'tail_prot')})

# --- per-variant summary (equilibrated portion) ---
summary = {}
print(f"\n{'Variant':<20s} {'EC50(pM)':>8s} | {'Tail-Prot mean±SD*':>18s} {'med':>5s} {'IQR':>5s} {'n_eff':>6s} {'CV':>7s} | "
      f"{'CA RMSD mean±SD*':>17s} {'med':>5s} {'IQR':>5s} {'n_eff':>6s} {'CV':>7s}")
print("  *mean±SD over replica means; med/IQR pooled over equilibrated frames; n_eff = sum of per-replica n_eff")
for v in VARIANTS:
    reps = sorted(results[v])
    if not reps: continue
    summary[v] = {}
    row = f"{LABELS[v]:<20s} {LAU_EC50[v]:>8.0f} | "
    for metric in ('tail_prot', 'ca_rmsd'):
        eq = {r: results[v][r][metric][results[v][r]['half']:] for r in reps}
        rep_means = [eq[r].mean() for r in reps]
        pooled = np.concatenate([eq[r] for r in reps])
        n_eff = sum(stats.summarize(eq[r])['n_eff'] for r in reps)
        summary[v][metric] = {
            'rep_means': [float(m) for m in rep_means],
            'mean': float(np.mean(rep_means)), 'sd': float(np.std(rep_means, ddof=1)) if len(reps) > 1 else 0.0,
            'median': float(np.median(pooled)),
            'iqr': float(np.percentile(pooled, 75) - np.percentile(pooled, 25)),
            'n_eff': float(n_eff), 'cv': stats.replica_cv(rep_means),
        }
        s = summary[v][metric]
        row += f"{s['mean']:>8.2f}±{s['sd']:<8.2f} {s['median']:>5.2f} {s['iqr']:>5.2f} {s['n_eff']:>6.0f} {s['cv']:>6.1%} | "
    print(row)

# --- significance tests ---
def pooled_eq(v, metric):
    """Concatenated equilibrated per-frame series; replicas truncated to common length."""
    reps = sorted(results[v])
    arrs = [results[v][r][metric][results[v][r]['half']:] for r in reps]
    L = min(len(a) for a in arrs)
    return np.concatenate([a[:L] for a in arrs])

tests = {'primary': [], 'secondary': [], 'control_repmeans': []}
print("\n=== PRIMARY: gglu_2oeg vs others, Tail-Prot, correlated t-test (concatenated per-frame) ===")
for other in VARIANTS:
    if other == "gglu_2oeg": continue
    x, y = pooled_eq("gglu_2oeg", 'tail_prot'), pooled_eq(other, 'tail_prot')
    L = min(len(x), len(y)); x, y = x[:L], y[:L]
    res = stats.correlated_t_test(x, y)
    res['other'] = other
    res['delta'] = float(x.mean() - y.mean())
    # control: unpaired t on per-replica means
    from scipy import stats as st
    m2 = summary['gglu_2oeg']['tail_prot']['rep_means']
    mo = summary[other]['tail_prot']['rep_means']
    tt = st.ttest_ind(m2, mo, equal_var=False)
    res['ctrl_t'], res['ctrl_df'], res['ctrl_p'] = float(tt.statistic), float(tt.df), float(tt.pvalue)
    tests['primary'].append(res)
    print(f"  2oeg vs {other:<11s} Δ={res['delta']:+.2f} Å  t={res['t']:.2f} df={res['df']:.0f} "
          f"p={res['p']:.3g} n_eff={res['n_eff']:.0f} | rep-mean Welch t={res['ctrl_t']:.2f} p={res['ctrl_p']:.3g}")

print("\n=== SECONDARY: gglu_3oeg vs others, CA RMSD, correlated t-test ===")
for other in VARIANTS:
    if other == "gglu_3oeg": continue
    x, y = pooled_eq("gglu_3oeg", 'ca_rmsd'), pooled_eq(other, 'ca_rmsd')
    L = min(len(x), len(y)); x, y = x[:L], y[:L]
    res = stats.correlated_t_test(x, y)
    res['other'] = other
    res['delta'] = float(x.mean() - y.mean())
    m3 = summary['gglu_3oeg']['ca_rmsd']['rep_means']
    mo = summary[other]['ca_rmsd']['rep_means']
    tt = st.ttest_ind(m3, mo, equal_var=False)
    res['ctrl_t'], res['ctrl_df'], res['ctrl_p'] = float(tt.statistic), float(tt.df), float(tt.pvalue)
    tests['secondary'].append(res)
    print(f"  3oeg vs {other:<11s} Δ={res['delta']:+.2f} Å  t={res['t']:.2f} df={res['df']:.0f} "
          f"p={res['p']:.3g} n_eff={res['n_eff']:.0f} | rep-mean Welch t={res['ctrl_t']:.2f} p={res['ctrl_p']:.3g}")

with open(f"{OUT}/test_results.json", "w") as f:
    json.dump({'summary': summary, 'tests': tests, 'skipped': skipped}, f, indent=1)
print(f"\nSaved {OUT}/test_results.json")

# --- faceted Tail-Prot distribution plot: rows = replica, cols = variant ---
fig, axes = plt.subplots(3, 4, figsize=(16, 9), sharex=True)
all_t = np.concatenate([results[v][r]['tail_prot'][results[v][r]['half']:]
                        for v in VARIANTS for r in results[v]])
bins = np.linspace(all_t.min(), all_t.max(), 60)
for j, v in enumerate(VARIANTS):
    for i, rep in enumerate([1, 2, 3]):
        ax = axes[i, j]
        if rep in results[v]:
            r = results[v][rep]
            d = r['tail_prot'][r['half']:]
            ax.hist(d, bins=bins, density=True, color=COLORS[v], alpha=0.75)
            ax.axvline(d.mean(), color='k', lw=1, ls='--')
            ax.text(0.95, 0.9, f"μ={d.mean():.2f}", transform=ax.transAxes,
                    ha='right', va='top', fontsize=8)
        else:
            ax.text(0.5, 0.5, "skipped", transform=ax.transAxes, ha='center')
        if i == 0: ax.set_title(LABELS[v], fontsize=10)
        if j == 0: ax.set_ylabel(f"rep{rep}\ndensity", fontsize=9)
        if i == 2: ax.set_xlabel("Tail-Prot min distance (Å)")
fig.suptitle("exp-D: C18 tail → protein min distance, equilibrated last 50 ns (dashed = replica mean)")
plt.tight_layout()
plt.savefig(f"{OUT}/tail_prot_compare.png", dpi=150)
print(f"Saved {OUT}/tail_prot_compare.png")
