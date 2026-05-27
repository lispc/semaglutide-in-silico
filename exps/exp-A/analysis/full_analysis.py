#!/usr/bin/env python3
"""Full 200 ns analysis: RMSD, catalytic distance, contacts, S1 pocket geometry, WT vs Aib8."""
import mdtraj as md, numpy as np, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

REPO = "/home/scroll/personal/semaglutide-in-silico"; OUT = f"{REPO}/exps/exp-A/analysis"
os.makedirs(OUT, exist_ok=True)

results = {}
for label in ["wt", "aib8"]:
    prmtop = f"{REPO}/exps/exp-A/tleap/wt.prmtop" if label == "wt" else f"{REPO}/exps/exp-A/md/aib8/aib8_modified.prmtop"
    dcd = f"{REPO}/exps/exp-A/md/{label}/{label}_traj.dcd"
    if not os.path.exists(dcd): continue

    print(f"Loading {label}...")
    t = md.load(dcd, top=prmtop)
    half = t.n_frames // 2
    t = t[half:]  # last 100 ns
    dt = 0.1  # ns per frame
    print(f"  {t.n_frames} frames ({t.n_frames*dt:.0f} ns)")

    # Selections
    dpp4_bb = t.topology.select("resid >= 0 and resid <= 727 and backbone")
    pep_bb  = t.topology.select("resid >= 728 and resid <= 758 and backbone")
    dpp4_heavy = t.topology.select("resid >= 0 and resid <= 727 and not element H")
    pep_heavy  = t.topology.select("resid >= 728 and resid <= 758 and not element H")
    # Peptide res2: Aib8's prmtop still labels CB2 as HA (ParmEd rename not persisted to prmtop file)
    pep2_c   = t.topology.select("resid == 729 and name C")
    pep2_cb  = t.topology.select("resid == 729 and name CB")
    pep2_cb2 = t.topology.select("resid == 729 and (name CB2 or name HA)")  # HA is actually CB2
    # Catalytic Ser
    ser_ogs = t.topology.select("resid >= 0 and resid <= 727 and resname SER and name OG")
    pn = t.topology.select("resid >= 728 and resid <= 733 and name CA")
    d = np.linalg.norm(t.xyz[0,ser_ogs]-t.xyz[0,pn].mean(0), axis=1)
    ser_og = ser_ogs[np.argmin(d)]
    # S1 pocket (heavy atoms)
    s1_heavy = t.topology.select(
        "(resid 547 or resid 629 or resid 631 or resid 656 or resid 662 or resid 666 or resid 710) and not element H")

    # Align
    t.superpose(t[0], atom_indices=dpp4_bb)

    # Metrics
    dpp4_rmsd = md.rmsd(t, t[0], atom_indices=dpp4_bb) * 10
    pep_rmsd  = md.rmsd(t, t[0], atom_indices=pep_bb) * 10
    cat_dist  = np.linalg.norm(t.xyz[:,ser_og]-t.xyz[:,pep2_c[0]], axis=1) * 10

    # Contacts between peptide heavy atoms and DPP-4 heavy atoms
    pep_xyz = t.xyz[:, pep_heavy]
    dpp4_xyz = t.xyz[:, dpp4_heavy]
    contacts_all = np.array([
        (np.linalg.norm(pep_xyz[i,:,np.newaxis]-dpp4_xyz[i,np.newaxis,:], axis=-1)*10 < 5.0).sum()
        for i in range(0, t.n_frames, 5)  # every 5th frame
    ])

    # S1 pocket: distances from CB/CB2 to pocket CA atoms
    s1_res = {547:'TYR', 629:'TRP', 631:'TYR', 656:'VAL', 662:'TYR', 666:'TYR', 710:'ASN'}
    s1_dists_cb = {}
    s1_dists_cb2 = {}
    for rid in s1_res:
        ca = t.topology.select(f"resid == {rid} and name CA")
        if len(ca) > 0 and len(pep2_cb) > 0:
            s1_dists_cb[rid] = np.linalg.norm(t.xyz[:,pep2_cb[0]]-t.xyz[:,ca[0]], axis=1)*10
        if len(pep2_cb2) > 0 and len(ca) > 0:
            s1_dists_cb2[rid] = np.linalg.norm(t.xyz[:,pep2_cb2[0]]-t.xyz[:,ca[0]], axis=1)*10

    results[label] = {
        'dpp4_rmsd': dpp4_rmsd, 'pep_rmsd': pep_rmsd, 'cat_dist': cat_dist,
        'contacts_all': contacts_all, 's1_dists_cb': s1_dists_cb, 's1_dists_cb2': s1_dists_cb2,
        'n_frames': t.n_frames, 'dt': dt, 'has_cb2': len(pep2_cb2) > 0
    }

# ====== SUMMARY ======
print(f"\n{'='*70}")
print(f"  FULL ANALYSIS: Last 100 ns of 200 ns trajectories")
print(f"{'='*70}")
print(f"  {'Metric':<35s} {'WT':>14s} {'Aib8':>14s} {'Δ':>10s}")
print(f"  {'-'*70}")

for metric, fmt, name in [
    ('dpp4_rmsd', '.1f', 'DPP-4 backbone RMSD (Å)'),
    ('pep_rmsd', '.1f', 'Peptide backbone RMSD (Å)'),
    ('cat_dist', '.1f', 'Catalytic distance (Å)'),
    ('contacts_all', '.0f', 'Peptide-DPP4 contacts (<5Å)'),
]:
    wm, ws = results['wt'][metric].mean(), results['wt'][metric].std()
    am, as_ = results['aib8'][metric].mean(), results['aib8'][metric].std()
    delta = am - wm
    print(f"  {name:<35s} {wm:{fmt}}±{ws:{fmt}} {am:{fmt}}±{as_:{fmt}} {delta:+{fmt}}")

# S1 pocket distances
print(f"  {'-'*70}")
print(f"  S1 Pocket Distances (CB/CB2 → residue CA, Å)")
for rid in [547, 629, 631, 656, 662, 666, 710]:
    rname = {547:'TYR',629:'TRP',631:'TYR',656:'VAL',662:'TYR',666:'TYR',710:'ASN'}[rid]
    w_cb = results['wt']['s1_dists_cb'].get(rid)
    a_cb = results['aib8']['s1_dists_cb'].get(rid)
    a_cb2 = results['aib8']['s1_dists_cb2'].get(rid)
    if w_cb is not None and a_cb is not None:
        print(f"  CB  → {rname}{rid:<4}           {w_cb.mean():>8.1f} ± {w_cb.std():.1f}  {a_cb.mean():>8.1f} ± {a_cb.std():.1f}")
    if a_cb2 is not None:
        print(f"  CB2 → {rname}{rid:<4}           {'--':>8s}          {a_cb2.mean():>8.1f} ± {a_cb2.std():.1f}")

# ====== PLOTS ======
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
time_wt = np.arange(results['wt']['n_frames']) * 0.1
time_aib = np.arange(results['aib8']['n_frames']) * 0.1

# RMSD comparison
ax = axes[0,0]
ax.plot(time_wt, results['wt']['dpp4_rmsd'], alpha=0.5, label='WT')
ax.plot(time_aib, results['aib8']['dpp4_rmsd'], alpha=0.5, label='Aib8')
ax.set_ylabel('Å'); ax.set_title(f"DPP-4 RMSD"); ax.legend()

ax = axes[0,1]
ax.plot(time_wt, results['wt']['pep_rmsd'], alpha=0.5, label='WT')
ax.plot(time_aib, results['aib8']['pep_rmsd'], alpha=0.5, label='Aib8')
ax.set_ylabel('Å'); ax.set_title("Peptide backbone RMSD"); ax.legend()

# Catalytic distance
ax = axes[0,2]
ax.plot(time_wt, results['wt']['cat_dist'], alpha=0.5, label=f'WT (μ={results["wt"]["cat_dist"].mean():.1f}Å)')
ax.plot(time_aib, results['aib8']['cat_dist'], alpha=0.5, label=f'Aib8 (μ={results["aib8"]["cat_dist"].mean():.1f}Å)')
ax.axhline(y=3.5, color='gray', linestyle='--', alpha=0.5)
ax.set_ylabel('Å'); ax.set_title("Ser630 OG → Ala8/Aib8 C"); ax.legend()

# Contacts
ax = axes[1,0]
ct_wt = results['wt']['contacts_all']
ct_aib = results['aib8']['contacts_all']
fw = np.arange(len(ct_wt)) * 0.1 * 5  # every 5th frame * 0.1 ns
fa = np.arange(len(ct_aib)) * 0.1 * 5
ax.plot(fw, ct_wt, alpha=0.5, label=f'WT (μ={ct_wt.mean():.0f})')
ax.plot(fa, ct_aib, alpha=0.5, label=f'Aib8 (μ={ct_aib.mean():.0f})')
ax.set_ylabel('Contacts'); ax.set_title("Peptide-DPP4 heavy-atom contacts <5Å"); ax.legend()

# Catalytic distance histogram
ax = axes[1,1]
ax.hist(results['wt']['cat_dist'], bins=50, alpha=0.5, density=True, label='WT')
ax.hist(results['aib8']['cat_dist'], bins=50, alpha=0.5, density=True, label='Aib8')
ax.axvline(x=3.5, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Å'); ax.set_title("Catalytic distance distribution"); ax.legend()

# S1 pocket distances
ax = axes[1,2]
rids = [547, 629, 631, 656, 662, 666, 710]
rnames = ['Y547','W629','Y631','V656','Y662','Y666','N710']
x = np.arange(len(rids))
w_vals = [results['wt']['s1_dists_cb'].get(r, np.array([0])).mean() for r in rids]
a_vals = [results['aib8']['s1_dists_cb'].get(r, np.array([0])).mean() for r in rids]
a2_vals = [results['aib8']['s1_dists_cb2'].get(r, np.array([np.nan])).mean() for r in rids]
w_err  = [results['wt']['s1_dists_cb'].get(r, np.array([0])).std() for r in rids]
a_err  = [results['aib8']['s1_dists_cb'].get(r, np.array([0])).std() for r in rids]
w = 0.25
ax.bar(x - w, w_vals, w, yerr=w_err, label='WT CB', alpha=0.7)
ax.bar(x, a_vals, w, yerr=a_err, label='Aib8 CB', alpha=0.7)
ax.bar(x + w, a2_vals, w, label='Aib8 CB2', alpha=0.7)
ax.set_xticks(x); ax.set_xticklabels(rnames, rotation=45)
ax.set_ylabel('Distance (Å)'); ax.set_title("S1 pocket: CB/CB2 → residue CA"); ax.legend(fontsize=8)

for ax in axes.flat: ax.set_xlabel("Time (ns)" if ax in [axes[0,0],axes[0,1],axes[0,2],axes[1,0]] else "")
plt.suptitle(f"exp-A: DPP-4 + GLP-1 WT vs Aib8 (last 100 ns)", fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{OUT}/full_analysis.png", dpi=120)
plt.close()
print(f"\n  Plot saved: {OUT}/full_analysis.png")
print(f"\nDone!")
