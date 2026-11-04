#!/usr/bin/env python3
"""Compare this tree's Stage-4 comparative result to Josei's independent
MPRAsnakeflow 0.6.2 + MPRAnalyze pass on the same HE016AH_0715 run.

Both passes use logFC = log2(Derived/Ancestral):
  logFC < 0 -> Ancestral-enriched      logFC > 0 -> Derived-enriched
Ours is empirically verified (see sign-convention check); Josei's is confirmed
by their own ancestral_hits.tsv / modern_hits.tsv split.

Join key: our locus 'seqNNNN' <- their 'libA_seqNNNN' / 'libB_seqNNNN'.
Verified 0 seq-ID collisions across libA/libB in MH_oligos.fa.
"""
import os, sys, math
import numpy as np
from scipy import stats

BASE = "/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo"
OURS = f"{BASE}/output/mpranalyze_comparative_rerun/results/adipo_comparative_rerun_final.txt"
JOSEI = "/wynton/group/ahituv/Josei_MPRA/downstream_mpranalyze/results/Derived_vs_Ancestral_all_globalFDR.tsv"
FDR = 0.05


def load_ours(path):
    d = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        # row.names written with col.names=NA -> first header field is empty
        cols = {c: i + 1 for i, c in enumerate(hdr[1:])} if hdr[0] == "" else {c: i for i, c in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            try:
                d[p[0]] = (float(p[cols["logFC"]]), float(p[cols["fdr"]]), float(p[cols["pval"]]))
            except (ValueError, KeyError, IndexError):
                continue
    return d


def load_josei(path):
    d = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(hdr):
                continue
            name = p[c["element"]]
            if "_seq" not in name:      # skips the one shell-escape-corrupted row
                continue
            key = name.split("_", 1)[1]          # libA_seq123 -> seq123
            try:
                d[key] = (float(p[c["logFC"]]), float(p[c["fdr"]]), float(p[c["pval"]]))
            except ValueError:
                continue
    return d


def main():
    if not os.path.exists(OURS):
        sys.exit(f"missing {OURS} - Stage 4 finalize has not run")
    ours, jos = load_ours(OURS), load_josei(JOSEI)

    o_sig = {k for k, v in ours.items() if v[1] < FDR}
    j_sig = {k for k, v in jos.items() if v[1] < FDR}

    def split(sig, src):
        a = sum(1 for k in sig if src[k][0] < 0)
        return a, len(sig) - a

    oa, od = split(o_sig, ours)
    ja, jd = split(j_sig, jos)

    print("=" * 78)
    print("MH adipocyte comparative (Ancestral vs Derived) -- two independent passes")
    print("=" * 78)
    print(f"\n{'':32}{'this tree':>14}{'Josei MPRAsnakeflow':>22}")
    print(f"{'loci tested':<32}{len(ours):>14,}{len(jos):>22,}")
    print(f"{'significant @ FDR<'+str(FDR):<32}{len(o_sig):>14,}{len(j_sig):>22,}")
    print(f"{'  Ancestral-enriched (logFC<0)':<32}{oa:>14,}{ja:>22,}")
    print(f"{'  Derived-enriched  (logFC>0)':<32}{od:>14,}{jd:>22,}")
    print(f"{'hit rate':<32}{100*len(o_sig)/max(len(ours),1):>13.1f}%{100*len(j_sig)/max(len(jos),1):>21.1f}%")

    shared = set(ours) & set(jos)
    print(f"\n--- shared loci: {len(shared):,} "
          f"({100*len(shared)/max(len(ours),1):.1f}% of ours, "
          f"{100*len(shared)/max(len(jos),1):.1f}% of theirs) ---")
    print(f"  ours-only: {len(set(ours)-set(jos)):,}   theirs-only: {len(set(jos)-set(ours)):,}")

    if not shared:
        sys.exit("\nNo shared loci -- check the name join.")

    ov_o = {k for k in shared if ours[k][1] < FDR}
    ov_j = {k for k in shared if jos[k][1] < FDR}
    both = ov_o & ov_j
    print(f"\n--- agreement on the {len(shared):,} shared loci ---")
    print(f"  significant in both      : {len(both):,}")
    print(f"  ours only                : {len(ov_o - ov_j):,}")
    print(f"  theirs only              : {len(ov_j - ov_o):,}")
    print(f"  neither                  : {len(shared) - len(ov_o | ov_j):,}")
    if ov_o:
        print(f"  recovery of our hits by theirs : {100*len(both)/len(ov_o):.1f}%")
    if ov_j:
        print(f"  recovery of their hits by ours : {100*len(both)/len(ov_j):.1f}%")

    # 2x2 Fisher on shared loci
    a = len(both); b = len(ov_o - ov_j); c = len(ov_j - ov_o)
    d = len(shared) - len(ov_o | ov_j)
    if min(a + b, a + c) > 0 and d > 0:
        orr, pv = stats.fisher_exact([[a, b], [c, d]])
        print(f"  Fisher 2x2 odds ratio {orr:.1f}  p={pv:.3g}")

    # direction concordance
    if both:
        conc = sum(1 for k in both if (ours[k][0] > 0) == (jos[k][0] > 0))
        print(f"\n  direction agreement among co-significant loci: "
              f"{conc}/{len(both)} = {100*conc/len(both):.1f}%")

    lo = np.array([ours[k][0] for k in shared])
    lj = np.array([jos[k][0] for k in shared])
    ok = np.isfinite(lo) & np.isfinite(lj)
    print(f"\n  logFC correlation over all shared loci (n={ok.sum():,}):")
    print(f"    Pearson  r = {stats.pearsonr(lo[ok], lj[ok])[0]:.3f}")
    print(f"    Spearman r = {stats.spearmanr(lo[ok], lj[ok])[0]:.3f}")
    print(f"    sign agreement = {100*((lo[ok]>0)==(lj[ok]>0)).mean():.1f}%")
    sl = np.array([k in ov_o for k in shared])[ok]
    if sl.sum() > 2:
        print(f"    restricted to our significant loci (n={sl.sum():,}): "
              f"Pearson r = {stats.pearsonr(lo[ok][sl], lj[ok][sl])[0]:.3f}, "
              f"sign agreement = {100*((lo[ok][sl]>0)==(lj[ok][sl]>0)).mean():.1f}%")

    print("\n--- our strongest hits, side by side ---")
    top = sorted(ov_o, key=lambda k: ours[k][1])[:15]
    print(f"  {'locus':<12}{'our logFC':>11}{'our FDR':>11}"
          f"{'their logFC':>13}{'their FDR':>12}  {'call'}")
    for k in top:
        ol, of_, _ = ours[k]; jl, jf, _ = jos[k]
        agree = "concordant" if (ol > 0) == (jl > 0) else "OPPOSITE SIGN"
        rep = "" if jf < FDR else " (n.s. in theirs)"
        print(f"  {k:<12}{ol:>11.4f}{of_:>11.2e}{jl:>13.4f}{jf:>12.2e}  {agree}{rep}")


if __name__ == "__main__":
    main()
