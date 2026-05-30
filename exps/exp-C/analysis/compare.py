#!/usr/bin/env python3
"""C18 monoacid vs diacid comparison: FA3 binding stability."""
import mdtraj as md, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

REPO = "/home/scroll/personal/semaglutide-in-silico"
OUT = f"{REPO}/exps/exp-C/analysis"
prmtop_m = f"{REPO}/exps/exp-C/tleap/c18_monoacid.prmtop"
prmtop_d = f"{REPO}/exps/exp-C/tleap/c18_diacid.prmtop"

data = {}
for label, prmtop, sys_name in [("monoacid", prmtop_m, "c18_monoacid"),
                                  ("diacid",  prmtop_d, "c18_diacid")]:
    data[label] = {}
    for rep in [1, 2, 3]:
        dcd = f"{REPO}/exps/exp-C/md/{sys_name}/rep{rep}/{sys_name}_traj.dcd"
        if not __import__('os').path.exists(dcd): continue
        t = md.load(dcd, top=prmtop)
        half = t.n_frames // 2; t_eq = t[half:]

        hsa_ca = t.topology.select("protein and name CA")
        fah_h = t.topology.select("resname FAH and not element H")
        o_distal = t.topology.select("resname FAH and (name O1D or name O2D)")
        o_proximal = t.topology.select("resname FAH and (name O1P or name O2P)")  # only diacid
        # Nearest ARG to distal carboxyl
        arg_near = t.topology.select("resSeq 482 and (name NH1 or name NH2)")
        # Also find second ARG near proximal carboxyl (for diacid)
        arg_sec = t.topology.select("resSeq 442 and (name NH1 or name NH2)")

        t.superpose(t[0], atom_indices=hsa_ca)

        hsa_rmsd = md.rmsd(t, t[0], atom_indices=hsa_ca) * 10
        fah_rmsd = md.rmsd(t, t[0], atom_indices=fah_h) * 10

        # Distances
        d_primary = []; d_secondary = []
        for frame in range(t.n_frames):
            if len(o_distal) > 0 and len(arg_near) > 0:
                od_c = t.xyz[frame, o_distal]; ar_c = t.xyz[frame, arg_near]
                d_primary.append(min(np.linalg.norm(o-n) for o in od_c for n in ar_c) * 10)
            if len(o_proximal) > 0 and len(arg_sec) > 0:
                op_c = t.xyz[frame, o_proximal]; as_c = t.xyz[frame, arg_sec]
                d_secondary.append(min(np.linalg.norm(o-n) for o in op_c for n in as_c) * 10)

        d_primary = np.array(d_primary); d_secondary = np.array(d_secondary)

        data[label][rep] = {
            'hsa_rmsd': hsa_rmsd, 'fah_rmsd': fah_rmsd,
            'd_primary': d_primary, 'd_secondary': d_secondary
        }

        eq_h = len(hsa_rmsd) // 2
        print(f"{label} rep{rep}: HSA={hsa_rmsd[eq_h:].mean():.1f}Å, "
              f"FA_RMSD={fah_rmsd[eq_h:].mean():.1f}Å, "
              f"Odist→ARG482={d_primary[eq_h:].mean():.1f}Å", end="")
        if len(d_secondary) > 0:
            print(f", Oprox→ARG442={d_secondary[eq_h:].mean():.1f}Å")
        else:
            print()

# === Summary table ===
print(f"\n{'='*75}")
print(f"  {'System':<12s} {'Rep':<5s} {'HSA RMSD':<12s} {'FA RMSD':<12s} {'1° anchor':<12s} {'2° anchor':<12s}")
print(f"  {'-'*73}")
for label in ["monoacid", "diacid"]:
    for rep in [1, 2, 3]:
        if rep not in data[label]: continue
        r = data[label][rep]
        h = len(r['hsa_rmsd']) // 2
        hsa = f"{r['hsa_rmsd'][h:].mean():.1f}±{r['hsa_rmsd'][h:].std():.1f}"
        fa = f"{r['fah_rmsd'][h:].mean():.1f}±{r['fah_rmsd'][h:].std():.1f}"
        p1 = f"{r['d_primary'][h:].mean():.1f}±{r['d_primary'][h:].std():.1f}"
        p2 = f"{r['d_secondary'][h:].mean():.1f}±{r['d_secondary'][h:].std():.1f}" if len(r['d_secondary'])>0 else "N/A"
        print(f"  {label:<12s} {rep:<5d} {hsa:<12s} {fa:<12s} {p1:<12s} {p2:<12s}")

# === Plots ===
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
colors = {'monoacid': 'blue', 'diacid': 'red'}
for col, metric_key in enumerate(['hsa_rmsd', 'fah_rmsd', 'd_primary']):
    for row, label in enumerate(['monoacid', 'diacid']):
        ax = axes[row, col]
        for rep in [1,2,3]:
            if rep not in data[label]: continue
            r = data[label][rep]
            vals = r[metric_key]
            t_ns = np.arange(len(vals)) * 0.1
            ax.plot(t_ns, vals, alpha=0.5, label=f'Rep {rep}')
        half_x = len(vals) * 0.05
        ax.axvline(x=half_x, color='gray', linestyle='--', alpha=0.3)
        ax.set_xlabel("Time (ns)")
        titles = {'hsa_rmsd': f'{label}: HSA RMSD (Å)', 'fah_rmsd': f'{label}: FA RMSD (Å)',
                   'd_primary': f'{label}: COO⁻ → ARG482 (Å)'}
        ax.set_title(titles.get(metric_key, metric_key))
        if col == 0: ax.legend(fontsize=7)

plt.suptitle("exp-C: C18 Monoacid vs Diacid — FA3 Binding Stability (100 ns ×3)", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{OUT}/compare.png", dpi=120)
plt.close()

# Diacid specific: secondary anchor distance plot
if 'diacid' in data and any(len(data['diacid'][r].get('d_secondary',[])) > 0 for r in [1,2,3]):
    fig, ax = plt.subplots(figsize=(10, 5))
    for rep in [1,2,3]:
        if rep not in data['diacid']: continue
        r = data['diacid'][rep]
        if len(r['d_secondary']) > 0:
            ax.plot(np.arange(len(r['d_secondary']))*0.1, r['d_secondary'], alpha=0.7, label=f'Rep {rep}')
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("Distance (Å)")
    ax.set_title("Diacid: Proximal COO⁻ → ARG442 distance")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"{OUT}/diacid_secondary.png", dpi=120)
    plt.close()

print(f"\nPlots: {OUT}/compare.png, {OUT}/diacid_secondary.png")
print("Done!")
