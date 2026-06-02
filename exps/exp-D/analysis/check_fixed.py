#!/usr/bin/env python3
"""Quick analysis of fixed MD: NZ-C bond, RMSD, tail distance."""
import mdtraj as md, numpy as np, os

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_D = f"{REPO}/exps/exp-D"

VARIANTS = ["no_linker", "gglu", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]
LABELS = {"no_linker": "No linker", "gglu": "gGlu only",
          "gglu_1oeg": "gGlu-1xOEG", "gglu_2oeg": "gGlu-2xOEG (Sema)", "gglu_3oeg": "gGlu-3xOEG"}
LAU_EC50 = {"no_linker": 269, "gglu": 9.9, "gglu_1oeg": 4.8, "gglu_2oeg": 6.2, "gglu_3oeg": 27.7}

results = {}
for vname in VARIANTS:
    prmtop = f"{EXP_D}/tleap/{vname}_fixed.prmtop"
    results[vname] = {}
    for rep in [1, 2, 3]:
        dcd = f"{EXP_D}/md/{vname}_fixed/rep{rep}/{vname}_fixed_traj.dcd"
        if not os.path.exists(dcd): continue
        chunks = [c for c in md.iterload(dcd, top=prmtop, chunk=100)]
        if not chunks: continue
        t = chunks[0] if len(chunks) == 1 else md.join(chunks)
        half = t.n_frames // 2

        prot_ca = t.topology.select("protein and name CA")
        lnk_c = t.topology.select("resname LNK and element C")
        nz = t.topology.select("resname LYS and name NZ")
        lnk_c0 = t.topology.select("resname LNK and name C")

        t.superpose(t[0], atom_indices=prot_ca)
        ca_rmsd = md.rmsd(t, t[0], atom_indices=prot_ca) * 10

        if len(nz) > 0 and len(lnk_c0) > 0:
            nz_c = np.linalg.norm(t.xyz[:, nz[0]] - t.xyz[:, lnk_c0[0]], axis=1) * 10
        else:
            nz_c = np.zeros(t.n_frames)

        tail_dist = np.zeros(t.n_frames)
        if len(lnk_c) > 2 and len(prot_ca) > 0:
            tail_c = lnk_c[2:]
            for f in range(t.n_frames):
                d = np.linalg.norm(t.xyz[f, tail_c][:, None] - t.xyz[f, prot_ca][None, :], axis=2)
                tail_dist[f] = d.min() * 10

        results[vname][rep] = {
            'ca_rmsd': ca_rmsd, 'nz_c': nz_c, 'tail_dist': tail_dist,
            'nframes': t.n_frames, 'half': half
        }
        h = half
        print(f"{vname} rep{rep}: {t.n_frames}frames  CA={ca_rmsd[h:].mean():.1f}+-{ca_rmsd[h:].std():.1f}A  "
              f"NZ-C={nz_c[h:].mean():.2f}+-{nz_c[h:].std():.2f}A  Tail={tail_dist[h:].mean():.1f}+-{tail_dist[h:].std():.1f}A")

print(f"\n{'='*90}")
print(f"  {'Variant':<20s} {'EC50':>6s} {'CA RMSD':>12s} {'NZ-C':>10s} {'Tail-Prot':>12s} {'Stability'}")
print(f"  {'-'*88}")
for vname in VARIANTS:
    reps = [r for r in results[vname].values()]
    if not reps: continue
    ca_avg = np.mean([r['ca_rmsd'][r['half']:].mean() for r in reps])
    ca_std = np.std([r['ca_rmsd'][r['half']:].mean() for r in reps])
    nz_avg = np.mean([r['nz_c'][r['half']:].mean() for r in reps])
    nz_std = np.std([r['nz_c'][r['half']:].mean() for r in reps])
    tail_avg = np.mean([r['tail_dist'][r['half']:].mean() for r in reps])
    tail_std = np.std([r['tail_dist'][r['half']:].mean() for r in reps])
    ec50 = LAU_EC50[vname]
    nz_ok = "OK" if 1.2 < nz_avg < 1.8 else "BROKEN"
    print(f"  {LABELS[vname]:<20s} {ec50:>6d}  {ca_avg:>5.1f}+-{ca_std:<4.1f}   {nz_avg:>5.2f}+-{nz_std:.2f}  {tail_avg:>5.1f}+-{tail_std:<4.1f}   {nz_ok}")
