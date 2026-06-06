#!/usr/bin/env python3
"""Compare 4 linker variants: ECD stability, C18-ECD distance, linker dynamics.

STATISTICAL LIMITATIONS (2026-06-06):
- Current script reports mean±std on autocorrelated MD frames. This overstates precision.
- Added median+IQR and effective sample size (n_eff) via common.lib.stats.
- Replica CV is now reported, but n=3 is still marginal for between-replica inference.
- No correlated t-test yet — significance claims must wait for n_eff-corrected testing.
- dt_ns is now derived from DCD timestamps instead of hard-coded 0.1 ns.
"""
import mdtraj as md, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os, sys

# Add project lib to path
REPO_ROOT = "/home/scroll/personal/semaglutide-in-silico"
sys.path.insert(0, f"{REPO_ROOT}/common/lib")
import stats

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_D = f"{REPO}/exps/exp-D"
TLEAP = f"{EXP_D}/tleap"
MD_DIR = f"{EXP_D}/md"
OUT = f"{EXP_D}/analysis"
os.makedirs(OUT, exist_ok=True)

VARIANTS = ["no_linker", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]
LABELS = {"no_linker": "No linker (19)", "gglu_1oeg": "γGlu-1×OEG (21)",
          "gglu_2oeg": "γGlu-2×OEG (Sema)", "gglu_3oeg": "γGlu-3×OEG (24)"}
COLORS = {"no_linker": "#E53935", "gglu_1oeg": "#FB8C00",
          "gglu_2oeg": "#43A047", "gglu_3oeg": "#1E88E5"}
LAU_EC50 = {"no_linker": 269, "gglu_1oeg": 4.8, "gglu_2oeg": 6.2, "gglu_3oeg": 27.7}

results = {}

for vname in VARIANTS:
    prmtop = f"{TLEAP}/{vname}.prmtop"
    results[vname] = {}
    for rep in [1, 2, 3]:
        dcd = f"{MD_DIR}/{vname}/rep{rep}/{vname}_traj.dcd"
        if not os.path.exists(dcd):
            print(f"  MISSING: {dcd}")
            continue

        print(f"Loading {vname} rep{rep}...")
        # Use iterload for incomplete DCD (MD just finished writing)
        chunks = []
        for chunk in md.iterload(dcd, top=prmtop, chunk=100):
            chunks.append(chunk)
        if not chunks:
            print(f"  EMPTY: {dcd}")
            continue
        t = chunks[0] if len(chunks) == 1 else md.join(chunks)
        print(f"  {t.n_frames} frames, {t.time[-1]:.0f} ps")
        # Use last 50% as equilibrated
        half = t.n_frames // 2
        t_eq = t[half:]

        # Atom selections
        hsa_ca = t.topology.select("protein and name CA")  # all protein CA (ECD + peptide)
        lnk_carbons = t.topology.select("resname LNK and element C")
        lnk_tail = lnk_carbons[2:] if len(lnk_carbons) > 2 else lnk_carbons

        # Align to protein CA
        if len(hsa_ca) == 0:
            print(f"  WARNING: no protein CA found, skipping")
            continue
        t.superpose(t[0], atom_indices=hsa_ca)

        # 1. Protein CA RMSD
        ca_rmsd = md.rmsd(t, t[0], atom_indices=hsa_ca) * 10

        # 2. C18 tail COM → protein min distance
        prot_xyz = t.xyz[:, hsa_ca]
        tail_xyz = t.xyz[:, lnk_tail] if len(lnk_tail) > 0 else np.zeros((t.n_frames, 1, 3))
        tail_com = tail_xyz.mean(axis=1)
        min_dists = []
        for f in range(t.n_frames):
            if len(lnk_tail) > 0:
                dists = np.linalg.norm(tail_xyz[f][:, None] - prot_xyz[f][None, :], axis=2)
                min_dists.append(dists.min() * 10)
            else:
                min_dists.append(0)
        tail_prot_dist = np.array(min_dists)

        # 4. Linker end-to-end: amide C to distal COO⁻ C
        # First C in tail is the first CH2 after amide
        # Last C in LNK is the COO⁻ carbon
        lnk_c_atoms = t.topology.select("resname LNK and element C")
        if len(lnk_c_atoms) >= 2:
            ee_start = lnk_c_atoms[1]   # first CH2 after amide
            ee_end = lnk_c_atoms[-1]     # last C (COO⁻ carbon)
            ee_vec = t.xyz[:, ee_end] - t.xyz[:, ee_start]
            ee_dist = np.linalg.norm(ee_vec, axis=1) * 10
        else:
            ee_dist = np.zeros(t.n_frames)

        # 5. NZ-C distance
        nz_atom = t.topology.select("resname LYS and name NZ")
        lnk_c0 = t.topology.select("resname LNK and name C")
        if len(nz_atom) > 0 and len(lnk_c0) > 0:
            nz_c = np.linalg.norm(t.xyz[:, nz_atom[0]] - t.xyz[:, lnk_c0[0]], axis=1) * 10
        else:
            nz_c = np.zeros(t.n_frames)

        # Derive dt from DCD timestamps instead of hard-coding
        dt_ns = (t.time[1] - t.time[0]) / 1000.0 if t.n_frames > 1 else 0.1

        results[vname][rep] = {
            'ca_rmsd': ca_rmsd, 'tail_prot_dist': tail_prot_dist,
            'ee_dist': ee_dist, 'nz_c': nz_c, 'half': half, 'nframes': t.n_frames,
            'dt_ns': dt_ns,
            'time_ns': np.arange(t.n_frames) * dt_ns
        }

        # Stats (equilibrated only) — upgraded with median, IQR, n_eff
        h = half
        for arr, label in [(ca_rmsd[h:], 'CA RMSD'), (tail_prot_dist[h:], 'Tail-prot'), (nz_c[h:], 'NZ-C')]:
            s = stats.summarize(arr, name=label)
            print(f"  {stats.format_summary(s)}")
        print(f"  (mean±std for compatibility: CA={ca_rmsd[h:].mean():.2f}±{ca_rmsd[h:].std():.2f}, "
              f"Tail={tail_prot_dist[h:].mean():.2f}±{tail_prot_dist[h:].std():.2f}, "
              f"NZ-C={nz_c[h:].mean():.2f}±{nz_c[h:].std():.2f})")

# === Per-replica summary (with n_eff) ===
print(f"\n{'='*100}")
print(f"  Per-replica stats (median + IQR + n_eff from autocorrelation correction)")
print(f"  {'-'*98}")
for vname in VARIANTS:
    for rep in [1, 2, 3]:
        if rep not in results[vname]: continue
        r = results[vname][rep]
        h = r['half']
        for arr, label in [(r['ca_rmsd'][h:], 'CA_RMSD'), (r['tail_prot_dist'][h:], 'Tail'), (r['nz_c'][h:], 'NZ-C')]:
            s = stats.summarize(arr, name=f"{LABELS[vname][:12]}_r{rep}_{label}")
            print(f"    {stats.format_summary(s)}")

# === Legacy mean±std summary (for cross-reference) ===
print(f"\n{'='*85}")
print(f"  {'Variant':<20s} {'Rep':>4s} {'CA RMSD':>12s} {'Tail-Prot':>10s} {'E2E':>10s} {'NZ-C':>10s}")
print(f"  {'-'*83}")
for vname in VARIANTS:
    for rep in [1, 2, 3]:
        if rep not in results[vname]: continue
        r = results[vname][rep]
        h = r['half']
        ecd = f"{r['ca_rmsd'][h:].mean():.1f}±{r['ca_rmsd'][h:].std():.1f}"
        pep = f"{r['ca_rmsd'][h:].mean():.1f}±{r['ca_rmsd'][h:].std():.1f}"
        te = f"{r['tail_prot_dist'][h:].mean():.1f}±{r['tail_prot_dist'][h:].std():.1f}"
        ee = f"{r['ee_dist'][h:].mean():.1f}±{r['ee_dist'][h:].std():.1f}"
        nc = f"{r['nz_c'][h:].mean():.1f}±{r['nz_c'][h:].std():.1f}"
        print(f"  {LABELS[vname]:<20s} {rep:>4d} {ecd:>12s} {pep:>12s} {te:>10s} {ee:>10s} {nc:>10s}")

# === Replica-averaged summary (with CV) ===
print(f"\n{'='*85}")
print(f"  {'Variant':<20s} {'EC50':>6s} {'ECD RMSD':>12s} {'Pep RMSD':>12s} {'Tail-ECD':>10s} {'Rep-CV':>8s}")
print(f"  {'-'*83}")
for vname in VARIANTS:
    reps = [r for r in results[vname].values()]
    if not reps: continue
    ecd_means = [r['ca_rmsd'][r['half']:].mean() for r in reps]
    te_means = [r['tail_prot_dist'][r['half']:].mean() for r in reps]
    ecd_avg = np.mean(ecd_means)
    ecd_std = np.std(ecd_means)
    te_avg = np.mean(te_means)
    te_std = np.std(te_means)
    cv = stats.replica_cv(te_means)  # CV on the most discriminating metric
    ec50 = float(LAU_EC50[vname])
    print(f"  {LABELS[vname]:<20s} {ec50:>6.0f} {ecd_avg:>6.1f}±{ecd_std:<4.1f}   {te_avg:>6.1f}±{te_std:<4.1f}   {cv:>7.1%}")

# === Plots ===
fig, axes = plt.subplots(2, 2, figsize=(16, 10))

for vname in VARIANTS:
    c = COLORS[vname]
    for rep in [1, 2, 3]:
        if rep not in results[vname]: continue
        r = results[vname][rep]
        t_ns = r['time_ns']
        axes[0,0].plot(t_ns, r['ca_rmsd'], color=c, alpha=0.4, linewidth=0.5)
        axes[0,1].plot(t_ns, r['tail_prot_dist'], color=c, alpha=0.4, linewidth=0.5)
        axes[1,0].plot(t_ns, r['ee_dist'], color=c, alpha=0.4, linewidth=0.5)
        axes[1,1].plot(t_ns, r['nz_c'], color=c, alpha=0.4, linewidth=0.5)
    # Bold mean line
    all_ecd = np.array([results[vname][r]['ca_rmsd'] for r in [1,2,3] if r in results[vname]])
    all_tail = np.array([results[vname][r]['tail_prot_dist'] for r in [1,2,3] if r in results[vname]])
    if len(all_ecd) > 0:
        axes[0,0].plot(t_ns, all_ecd.mean(0), color=c, linewidth=2, label=LABELS[vname])
        axes[0,1].plot(t_ns, all_tail.mean(0), color=c, linewidth=2, label=LABELS[vname])

axes[0,0].set_ylabel('ECD CA RMSD (Å)'); axes[0,0].set_title('ECD Backbone Stability')
axes[0,1].set_ylabel('Distance (Å)'); axes[0,1].set_title('C18 Tail → ECD Min Distance')
axes[1,0].set_ylabel('Distance (Å)'); axes[1,0].set_title('Linker End-to-End Distance')
axes[1,1].set_ylabel('Distance (Å)'); axes[1,1].set_title('Lys26 NZ – Amide C Distance')
for ax in axes.flat: ax.legend(fontsize=8); ax.set_xlabel('Time (ns)')
plt.tight_layout()
plt.savefig(f"{OUT}/linker_compare.png", dpi=150)
print(f"\nSaved: {OUT}/linker_compare.png")

# === Bar chart comparison ===
fig2, axes2 = plt.subplots(1, 3, figsize=(14, 5))
vnames = [v for v in VARIANTS if results[v]]
x = np.arange(len(vnames))
ecd_means = [np.mean([results[v][r]['ca_rmsd'][results[v][r]['half']:].mean() for r in [1,2,3] if r in results[v]]) for v in vnames]
ecd_stds = [np.std([results[v][r]['ca_rmsd'][results[v][r]['half']:].mean() for r in [1,2,3] if r in results[v]]) for v in vnames]
tail_means = [np.mean([results[v][r]['tail_prot_dist'][results[v][r]['half']:].mean() for r in [1,2,3] if r in results[v]]) for v in vnames]
tail_stds = [np.std([results[v][r]['tail_prot_dist'][results[v][r]['half']:].mean() for r in [1,2,3] if r in results[v]]) for v in vnames]
ee_means = [np.mean([results[v][r]['ee_dist'][results[v][r]['half']:].mean() for r in [1,2,3] if r in results[v]]) for v in vnames]
ee_stds = [np.std([results[v][r]['ee_dist'][results[v][r]['half']:].mean() for r in [1,2,3] if r in results[v]]) for v in vnames]

colors_list = [COLORS[v] for v in vnames]
axes2[0].bar(x, ecd_means, yerr=ecd_stds, color=colors_list, capsize=5)
axes2[0].set_xticks(x); axes2[0].set_xticklabels([LABELS[v] for v in vnames], rotation=15, ha='right', fontsize=8)
axes2[0].set_ylabel('ECD RMSD (Å)'); axes2[0].set_title('ECD Stability')
axes2[1].bar(x, tail_means, yerr=tail_stds, color=colors_list, capsize=5)
axes2[1].set_xticks(x); axes2[1].set_xticklabels([LABELS[v] for v in vnames], rotation=15, ha='right', fontsize=8)
axes2[1].set_ylabel('Distance (Å)'); axes2[1].set_title('C18 Tail – ECD Distance')
axes2[2].bar(x, ee_means, yerr=ee_stds, color=colors_list, capsize=5)
axes2[2].set_xticks(x); axes2[2].set_xticklabels([LABELS[v] for v in vnames], rotation=15, ha='right', fontsize=8)
axes2[2].set_ylabel('Distance (Å)'); axes2[2].set_title('Linker End-to-End')
plt.tight_layout()
plt.savefig(f"{OUT}/linker_bars.png", dpi=150)
print(f"Saved: {OUT}/linker_bars.png")
