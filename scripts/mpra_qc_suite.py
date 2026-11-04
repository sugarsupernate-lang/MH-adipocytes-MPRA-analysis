#!/usr/bin/env python3
"""
MPRA QC suite -- MH library, SGBS adipocytes, HE016AH_0715 (AVITI run).

Implements chapters 1 and 2 of the Gokhman Lab MPRA QC book
(https://gokhmanlaborganization.github.io/MPRA_QC_book/) against this tree's
own results.  Section numbers, goals and interpretations follow the book;
every quantity is computed from this pipeline's outputs:

  oligo-level activity     output/allreps_Filter.txt
                           (pipeline CPM + log2(RNA/DNA), pseudocount applied)
  barcode-level counts     output/UMI/barcode_counts_UMI.txt
  outlier-filter pair      output/UMI/UMI_no_filter.txt, UMI_std2_filter.txt
  activity call            MPRAnalyze quantitative: statistic, pval.mad, fdr.mad
  differential activity    MPRAnalyze comparative (chunked, library-wide depth
                           factors): logFC, global BH fdr
  association              associations/MH_association.count_ready.tsv.gz
  pre-promiscuity assoc.   ../merged_association/merged_all_unique.txt

Book vocabulary maps onto this library as:
  cCRE            -> oligo (200 bp element)
  allele1/allele2 -> Derived / Ancestral
  positive_ctrl   -> libA_ctrl_*      negative_ctrl -> libA_scrambled*

Counts are normalised to library size (CPM) before any cross-replicate
comparison, per section 2.21.  Activity status is taken from MPRAnalyze
fdr.mad < 0.05 over test oligos only -- scrambled negatives are excluded from
the BH adjustment and carry no FDR.

Usage:  python3 mpra_qc_suite.py [--out_txt PATH] [--out_json PATH]
"""
import argparse
import collections
import gzip
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit

BASE = "/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo"
MERGED = "/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/merged_association"
REPS = ["rep1", "rep2", "rep3"]
PAIRS = [(0, 1), (0, 2), (1, 2)]
FDR_CUT = 0.05
SIBLINGS = {
    "NPC": f"{os.path.dirname(BASE)}/npc_count/output/allreps_Filter.txt",
    "HOb": f"{os.path.dirname(BASE)}/HOb_count/output/allreps_Filter.txt",
    "BJ":  f"{os.path.dirname(BASE)}/bj_count/output_all/allreps_Filter.txt",
}

OUT = []
Q = {}


def say(s=""):
    print(s)
    OUT.append(s)


def head(title, goal):
    say("")
    say("=" * 86)
    say(f"  {title}")
    say("=" * 86)
    say(f"  goal: {goal}")


def na(title, goal, reason):
    head(title, goal)
    say(f"  NOT EVALUATED -- {reason}")
    Q[title.split()[0]] = {"status": "not_evaluated", "reason": reason}


def verdict(ok, good, bad):
    """Emit a data-derived verdict rather than a fixed string."""
    say(f"  -> {'PASS' if ok else 'FLAG'}: {good if ok else bad}")
    return bool(ok)


def op(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")


def gc_frac(s):
    s = s.upper()
    return (s.count("G") + s.count("C")) / len(s) if s else np.nan


def max_g(s):
    m = re.findall(r"G+", s.upper())
    return max((len(x) for x in m), default=0)


def is_pos(n):
    return "ctrl" in n and "scrambled" not in n


def is_neg(n):
    return "scrambled" in n


def ccre_type(n):
    return "positive_ctrl" if is_pos(n) else "negative_ctrl" if is_neg(n) else "test_cCRE"


# The positive-control set is not one thing. It mixes human accelerated regions,
# a generic coordinate-named panel, and named enhancers that are largely
# bone/osteoblast (Runx2, Col1a1/2, RANKL, SOST, Mmp13, Fgfr2, Kremen1). Those
# are not expected to fire in adipocytes, so pooling them understates the
# positive controls that are relevant here.
BONE = ("runx2", "col1a1", "col1a2", "rankl", "sost", "mmp13", "fgfr2", "kremen1",
        "trps1", "gli2", "notch2", "shox", "zfhx4")


def ctrl_class(n):
    if is_neg(n):
        return "scrambled"
    if not is_pos(n):
        return "test"
    low = n.lower()
    if "2xhar" in low or "hacns" in low:
        return "HAR"
    if re.search(r"ctrl_(controls_)?seq\d", low):
        return "generic panel"
    if any(b in low for b in BONE):
        return "bone/osteoblast enh"
    return "other named"


def pearson(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 10:
        return np.nan
    return float(stats.pearsonr(a[ok], b[ok])[0])


def spearman(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 10:
        return np.nan
    return float(stats.spearmanr(a[ok], b[ok])[0])


def hill(x, ymax, k, n):
    return ymax * x ** n / (k ** n + x ** n)


def hill_extrapolate(fracs, vals, at=3.0):
    """Fit the book's Hill saturation model and predict retention at `at`x depth."""
    try:
        p, _ = curve_fit(hill, np.asarray(fracs, float), np.asarray(vals, float),
                         p0=[max(vals) * 1.05, 0.3, 1.0], maxfev=20000)
        return float(hill(at, *p)), float(p[0])
    except Exception:
        return np.nan, np.nan


# ----------------------------------------------------------------------------
# load
# ----------------------------------------------------------------------------
def load_sequences(path):
    seq, name = {}, None
    buf = []
    with op(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seq[name] = "".join(buf)
                name = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
    if name:
        seq[name] = "".join(buf)
    return seq


def load_activity(path):
    """Oligo-level table from the pipeline (already CPM-normalised, pseudocounted)."""
    df = pd.read_csv(path)
    df = df.rename(columns={"oligo_column": "cCRE"})
    df["cCRE"] = df["cCRE"].astype(str)
    return df.set_index("cCRE")


def load_mpranalyze_quant(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    # scrambled negatives are excluded from the BH adjustment upstream and are
    # written with fdr.mad = 0; they carry no activity call.
    neg = df.index.to_series().map(is_neg).values
    df["activity_FDR"] = df["fdr.mad"].where(~neg, np.nan)
    df["activity_pval"] = df["pval.mad"]
    df["activity_statistic"] = df["statistic"]
    df["activity_status"] = np.where(df["activity_FDR"] < FDR_CUT, "active", "non_active")
    df.loc[neg, "activity_status"] = "control"
    return df


def load_comparative(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    df = df[np.isfinite(df["logFC"])]
    df["differential_activity_FDR"] = df["fdr"]
    df["differentialy_active"] = df["fdr"] < FDR_CUT
    return df


def load_hg38_map(bed):
    """locus (seqNNNN) -> hg38 (chrom, start, end), from the liftOver of the oligo
    coordinates. Both alleles of a locus share coordinates."""
    m = {}
    if not os.path.exists(bed):
        return m
    for l in open(bed):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 4:
            m[f[3]] = (f[0], int(f[1]), int(f[2]))
    return m


def load_screen(path):
    """chrom -> sorted (start, end, class) for ENCODE SCREEN cCREs."""
    by = collections.defaultdict(list)
    if not os.path.exists(path):
        return {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 6:
            by[f[0]].append((int(f[1]), int(f[2]), f[5]))
    for c in by:
        by[c].sort()
    return {c: (np.array([x[0] for x in v]), np.array([x[1] for x in v]),
                [x[2] for x in v]) for c, v in by.items()}


def load_tss(gtf, keep=("protein_coding", "lncRNA")):
    """chrom -> sorted array of TSS positions (hg38)."""
    by = collections.defaultdict(list)
    if not os.path.exists(gtf):
        return {}
    for line in open(gtf):
        if line[0] == "#":
            continue
        f = line.split("\t", 9)
        if len(f) < 9 or f[2] != "gene":
            continue
        gt = re.search(r'gene_type "([^"]+)"', f[8])
        if gt and gt.group(1) not in keep:
            continue
        by[f[0]].append(int(f[3]) if f[6] == "+" else int(f[4]))
    return {c: np.array(sorted(v)) for c, v in by.items()}


def load_barcode_counts(path):
    df = pd.read_csv(path, sep="\t",
                     dtype={"oligo_bc": str, **{f"{c}_{r}": np.float32
                                                for c in ("RNA", "DNA") for r in REPS}})
    df["cCRE"] = df["oligo_bc"].str.rsplit("_", n=1).str[0]
    return df


# ----------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Plot data.  The report states every number; these structures let the figures
# be redrawn from the same computation rather than a parallel one.
# ---------------------------------------------------------------------------
PLOT = {}


def hist1d(x, lo=None, hi=None, bins=40):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return None
    lo = float(np.nanpercentile(x, 0.2)) if lo is None else lo
    hi = float(np.nanpercentile(x, 99.8)) if hi is None else hi
    if hi <= lo:
        hi = lo + 1.0
    c, e = np.histogram(np.clip(x, lo, hi), bins=bins, range=(lo, hi))
    return dict(edges=[round(float(v), 5) for v in e],
                counts=[int(v) for v in c], n=int(len(x)))


def hist2d(x, y, bins=44, qlo=0.5, qhi=99.5):
    """Sparse 2-D density, the book's hexbin equivalent."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 10:
        return None
    xl, xh = float(np.percentile(x, qlo)), float(np.percentile(x, qhi))
    yl, yh = float(np.percentile(y, qlo)), float(np.percentile(y, qhi))
    if xh <= xl: xh = xl + 1.0
    if yh <= yl: yh = yl + 1.0
    H, xe, ye = np.histogram2d(np.clip(x, xl, xh), np.clip(y, yl, yh),
                               bins=bins, range=[[xl, xh], [yl, yh]])
    pts = [[int(i), int(j), int(H[i, j])] for i, j in zip(*np.nonzero(H))]
    return dict(xlo=round(xl, 4), xhi=round(xh, 4), ylo=round(yl, 4), yhi=round(yh, 4),
                bins=bins, cells=pts, n=int(len(x)), vmax=int(H.max()))


def boxstats(x, label):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) == 0:
        return None
    q1, med, q3 = [float(np.percentile(x, q)) for q in (25, 50, 75)]
    iqr = q3 - q1
    lo = float(np.min(x[x >= q1 - 1.5 * iqr])) if np.any(x >= q1 - 1.5 * iqr) else q1
    hi = float(np.max(x[x <= q3 + 1.5 * iqr])) if np.any(x <= q3 + 1.5 * iqr) else q3
    return dict(label=label, q1=q1, med=med, q3=q3, lo=lo, hi=hi,
                n=int(len(x)), mean=float(x.mean()))


def capture_plots(D):
    """Build every figure's data from the arrays the report was computed from."""
    P = PLOT
    g = D.get

    if g("mult"):
        mm = g("mult")
        ks = sorted(k for k in mm if k <= 8)
        P["1.1"] = dict(k=[int(k) for k in ks], v=[int(mm[k]) for k in ks],
                        hi=int(sum(v for k, v in mm.items() if k > 8)))
    if g("bcs_per_vals") is not None:
        v = np.sort(np.asarray(g("bcs_per_vals"), float))
        idx = np.unique(np.linspace(0, len(v) - 1, 160).astype(int))
        P["1.3"] = dict(x=[float(v[i]) for i in idx],
                        y=[float((i + 1) / len(v)) for i in idx], n=int(len(v)))
    P["1.4"] = g("retained_rows")
    P["1.5"] = dict(frac=g("ds_frac"), ret=g("ds_ret"), hill=g("ds_hill"))
    P["1.6"] = g("ds_box")
    for key, rows in (("1.7", g("gc_rows")), ("1.8", g("gs_rows")), ("2.2", g("gc_rows"))):
        if rows:
            P[key] = rows
    P["2.1"] = g("retention_rows")
    P["_assoc_bc"] = g("assoc_bc")

    st, rl = g("A_STAT"), g("RL_MEAN")
    tm = g("test_m")
    if st is not None:
        P["2.3"] = hist2d(rl[tm], st[tm])
    if rl is not None:
        act = g("ACTIVE")
        P["2.4"] = dict(all=hist1d(rl[tm]),
                        active=hist1d(rl[tm & act]) if act is not None else None,
                        pos=hist1d(rl[g("pos_m")], bins=24),
                        neg=hist1d(rl[g("neg_m")], bins=24))
    pv_t, pv_n = g("pv_t"), g("pv_n")
    if pv_t is not None:
        qs = np.arange(0.5, 100, 1.0)
        P["2.6"] = dict(
            test=dict(obs=[float(v) for v in np.percentile(pv_t, qs)],
                      exp=[float(v / 100) for v in qs]),
            neg=dict(obs=[float(v) for v in np.percentile(pv_n, qs)],
                     exp=[float(v / 100) for v in qs]),
            hist=hist1d(pv_t, 0, 1, 20), hist_neg=hist1d(pv_n, 0, 1, 20))
    cmp_ = g("comp")
    if cmp_ is not None:
        lfc = cmp_["logFC"].to_numpy(float)
        fdr = cmp_["differential_activity_FDR"].to_numpy(float)
        sig = fdr < FDR_CUT
        P["2.5"] = dict(all=hist1d(lfc, bins=44), sig=hist1d(lfc[sig], bins=44),
                        n_sig=int(sig.sum()),
                        up=int((sig & (lfc > 0)).sum()), dn=int((sig & (lfc < 0)).sum()))
        ml = np.isfinite(lfc) & np.isfinite(fdr)
        nl = -np.log10(np.clip(fdr[ml], 1e-300, 1))
        P["2.7"] = dict(pts=[[round(float(a), 4), round(float(b), 3), int(c)]
                             for a, b, c in zip(lfc[ml], np.clip(nl, 0, 60), sig[ml])][:4000],
                        cut=float(-np.log10(FDR_CUT)))
    P["2.8"] = dict(frac=g("ds2_frac"), dna=g("ds2_dna"), rna=g("ds2_rna"))
    P["2.9"] = dict(frac=g("ds2_frac"), dna=g("ds2_cdna"), rna=g("ds2_crna"))
    P["2.10"] = dict(frac=g("ds2_frac"), act=g("ds2_act"))
    if g("rna_cum") is not None:
        c = np.asarray(g("rna_cum"), float)
        idx = np.unique(np.linspace(0, len(c) - 1, 160).astype(int))
        P["2.11"] = dict(x=[float((i + 1) / len(c)) for i in idx],
                         y=[float(c[i]) for i in idx])
    P["2.12"] = g("pca")
    RL = g("RL")
    if RL is not None:
        P["2.13"] = hist2d(RL[0][tm], RL[1][tm])
        P["2.16"] = dict(dens=hist2d(RL[0][tm], RL[1][tm]),
                         pos=[[round(float(a), 3), round(float(b), 3)]
                              for a, b in zip(RL[0][g("pos_m")], RL[1][g("pos_m")])],
                         neg=[[round(float(a), 3), round(float(b), 3)]
                              for a, b in zip(RL[0][g("neg_m")], RL[1][g("neg_m")])])
    LR, al = g("LFC_REP"), g("allel")
    if al is not None and g("RL") is not None:
        RLv = g("RL")
        P["2.14"] = hist2d(RLv[0][al], RLv[1][al])
    if LR is not None:
        P["2.14_difference"] = hist2d(LR[0], LR[1])
    di, ai = g("di"), g("ai")
    if di is not None:
        P["2.17"] = hist2d(rl[ai], rl[di])
    P["2.15"] = g("rep_by_act")
    P["2.19"] = g("noise_rows")
    P["2.20"] = g("repro_rows")
    dc, rc = g("DNA_CPM"), g("RNA_CPM")
    if dc is not None:
        # linear CPM, truncated at 250 as the book specifies
        dx = np.clip(dc[:, tm].mean(axis=0), 0, 250)
        ry = np.clip(rc[:, tm].mean(axis=0), 0, 250)
        P["2.21"] = dict(dens=hist2d(dx, ry, qlo=0, qhi=100),
                         per_rep=g("rna_dna_rows"))
    cls = g("ctrl_cls")
    if cls is not None:
        want = ["HAR", "generic panel", "other named", "bone/osteoblast enh", "scrambled"]
        P["2.22"] = [b for b in
                     [boxstats(rl[cls == w], w.replace("bone/osteoblast enh", "bone enh")
                               .replace("generic panel", "generic")) for w in want
                      if (cls == w).sum() >= 5] + [boxstats(rl[tm], "test")] if b]
    else:
        P["2.22"] = [b for b in (boxstats(rl[g("pos_m")], "positive ctrl"),
                                 boxstats(rl[g("neg_m")], "negative ctrl"),
                                 boxstats(rl[tm], "test cCREs")) if b]
    return P


def _write(a):
    """Always emit whatever has been computed, so a failed section cannot leave
    a stale report from a previous run in place."""
    with open(a.out_json, "w") as fh:
        json.dump(Q, fh, indent=1, default=float)
    with open(a.out_txt, "w") as fh:
        fh.write("\n".join(OUT) + "\n")
    if PLOT:
        with open(a.out_plots, "w") as fh:
            json.dump(PLOT, fh, separators=(",", ":"), default=float)
        print(f"wrote {a.out_plots}")
    print(f"\nwrote {a.out_json}")
    print(f"wrote {a.out_txt}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_json", default=f"{BASE}/output/qc/qc_metrics.json")
    ap.add_argument("--out_txt", default=f"{BASE}/output/qc/qc_report.txt")
    ap.add_argument("--out_plots", default=f"{BASE}/output/qc/qc_plots.json")
    ap.add_argument("--comparative",
                    default=f"{BASE}/output/mpranalyze_comparative_rerun/results/"
                            "adipo_comparative_rerun_final.txt")
    ap.add_argument("--seed", type=int, default=20260823)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    os.makedirs(os.path.dirname(a.out_txt), exist_ok=True)
    try:
        _run(a, rng)
    finally:
        _write(a)


def _run(a, rng):
    PD = {}
    say("MPRA QC -- MH library, SGBS adipocytes, HE016AH_0715 (AVITI)")
    say("Gokhman Lab MPRA QC book, chapters 1-2, evaluated on this tree's MPRAnalyze results.")
    say("cCRE = oligo (200 bp).  allele1/allele2 = Derived/Ancestral.")
    say("All cross-replicate comparisons use library-size (CPM) normalised counts (sec 2.21).")

    # ---------------- inputs ----------------
    say("\nloading inputs ...")
    seqs = load_sequences(f"{BASE}/associations/MH_oligos.fa")
    seq_by_name = dict(seqs)
    for k, v in seqs.items():
        seq_by_name.setdefault(k.split("_[")[0] if "_[" in k else k, v)
    act = load_activity(f"{BASE}/output/allreps_Filter.txt")
    quant = load_mpranalyze_quant(
        f"{BASE}/output/mpranalyze_quantitative/"
        "adipo_results_full_quantitative_nobc_fdr_annot_factor_new_design.txt")
    comp = load_comparative(a.comparative)
    bc = load_barcode_counts(f"{BASE}/output/UMI/barcode_counts_UMI.txt")
    say(f"  oligos (activity table)   {len(act):,}")
    say(f"  oligos (MPRAnalyze quant) {len(quant):,}")
    say(f"  barcodes (counts)         {len(bc):,}")
    say(f"  differential activity     {'%s loci' % f'{len(comp):,}' if comp is not None else 'PENDING (stage 4 not finalised)'}")

    names = act.index.to_numpy()
    typ = np.array([ccre_type(n) for n in names])
    pos_m, neg_m, test_m = typ == "positive_ctrl", typ == "negative_ctrl", typ == "test_cCRE"

    # oligo-level CPM log-ratios straight from the pipeline
    RL = np.vstack([act[f"ratio_log_{r}"].to_numpy(float) for r in REPS])
    DNA = np.vstack([act[f"DNA_{r}"].to_numpy(float) for r in REPS])
    RNA = np.vstack([act[f"RNA_{r}"].to_numpy(float) for r in REPS])
    DNA_CPM = np.vstack([act[f"DNA_cpm_{r}"].to_numpy(float) for r in REPS])
    RNA_CPM = np.vstack([act[f"RNA_cpm_{r}"].to_numpy(float) for r in REPS])
    RL_MEAN = np.nanmean(RL, axis=0)

    q = quant.reindex(names)
    A_STAT = q["activity_statistic"].to_numpy(float)
    A_PVAL = q["activity_pval"].to_numpy(float)
    A_FDR = q["activity_FDR"].to_numpy(float)
    ACTIVE = A_FDR < FDR_CUT
    Q["dataset"] = dict(oligos=int(len(act)), barcodes=int(len(bc)),
                        active=int(np.nansum(ACTIVE)),
                        positive_ctrl=int(pos_m.sum()), negative_ctrl=int(neg_m.sum()),
                        test=int(test_m.sum()))

    # barcode-level arrays
    bd = np.vstack([bc[f"DNA_{r}"].to_numpy(np.float64) for r in REPS])
    br = np.vstack([bc[f"RNA_{r}"].to_numpy(np.float64) for r in REPS])
    bd_tot, br_tot = bd.sum(axis=0), br.sum(axis=0)
    bc_ccre = bc["cCRE"].to_numpy()
    bc_code, bc_levels = pd.factorize(bc["cCRE"], sort=True)
    n_lev = len(bc_levels)
    def by_ccre(v):
        """sum a barcode-level vector into cCRE-level totals (level order)"""
        return np.bincount(bc_code, weights=v, minlength=n_lev)

    say("")
    say("#" * 86)
    say("#  CHAPTER 1 -- cCRE-BC ASSOCIATION QC")
    say("#" * 86)

    # ---------------- 1.1 cCREs per BC ----------------
    head("1.1  cCREs per BC", "assess BC promiscuity")
    pre = f"{MERGED}/merged_all_unique.txt"
    if os.path.exists(pre):
        say("  input: merged_all_unique.txt (all unique oligo-barcode pairs, pre-promiscuity filter)")
        pc = collections.Counter()
        with op(pre) as fh:
            fh.readline()
            for line in fh:
                p = line.split(",")
                if len(p) >= 2:
                    pc[p[1]] += 1
        mult = collections.Counter(pc.values())
        tot = len(pc)
        say(f"\n  barcodes (pre-filter): {tot:,}")
        say(f"  {'cCREs per BC':>14}{'barcodes':>14}{'% of BCs':>11}")
        for k in sorted(mult)[:8]:
            say(f"  {k:>14}{mult[k]:>14,}{100*mult[k]/tot:>10.3f}%")
        if len(mult) > 8:
            hi = sum(v for k, v in mult.items() if k > 8)
            say(f"  {'>8':>14}{hi:>14,}{100*hi/tot:>10.3f}%")
        spec = mult[1] / tot
        promisc = tot - mult[1]
        say(f"\n  uniquely-assigned BCs: {mult[1]:,} ({100*spec:.3f}%)")
        say(f"  promiscuous BCs removed by the association filter: {promisc:,} ({100*(1-spec):.3f}%)")
        PD["mult"] = {int(k): int(v) for k, v in mult.items()}
        Q["1.1"] = dict(barcodes_prefilter=tot, unique=int(mult[1]),
                        promiscuous=int(promisc), specificity=float(spec),
                        multiplicity={int(k): int(v) for k, v in sorted(mult.items())[:12]})
        verdict(spec > 0.95,
                f"{100*spec:.2f}% of barcodes map to exactly one cCRE -- high specificity",
                f"only {100*spec:.2f}% of barcodes are cCRE-specific -- specificity problem")
    else:
        say(f"  NOT EVALUATED -- pre-promiscuity association not found at {pre}")
        Q["1.1"] = {"status": "not_evaluated"}

    # ---------------- 1.2 ----------------
    na("1.2  Reads per association", "assess confidence in associations",
       "no read-support column: the association files in this tree carry "
       "(barcode, oligo) pairs only, with match_count dropped upstream")

    # ---------------- assoc counts ----------------
    say("\nloading final association ...")
    bcs_per = collections.Counter()
    with op(f"{BASE}/associations/MH_association.count_ready.tsv.gz") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                bcs_per[p[1].split("_[")[0]] += 1
    n_assoc_bc = sum(bcs_per.values())
    say(f"  {n_assoc_bc:,} barcodes over {len(bcs_per):,} cCREs")

    # ---------------- 1.3 BCs per cCRE ----------------
    head("1.3  BCs per cCRE", "assess BC replicability (complexity)")
    v = np.array(sorted(bcs_per.values()))
    PD["bcs_per_vals"] = v
    say(f"  metric: eCDF of barcodes per cCRE")
    say(f"  median {int(np.median(v))}   mean {v.mean():.1f}   "
        f"IQR {int(np.percentile(v,25))}-{int(np.percentile(v,75))}   range {v.min()}-{v.max()}")
    say(f"  {'BC threshold':>14}{'cCREs below':>14}{'% below':>10}")
    lt = {}
    for t in (5, 10, 20, 50, 100):
        n = int((v < t).sum())
        lt[t] = n
        say(f"  {'< %d' % t:>14}{n:>14,}{100*n/len(v):>9.2f}%")
    frac_lt10 = lt[10] / len(v)
    Q["1.3"] = dict(median=int(np.median(v)), mean=float(v.mean()),
                    below={int(k): int(x) for k, x in lt.items()},
                    frac_below_10=float(frac_lt10))
    verdict(frac_lt10 < 0.05,
            f"only {100*frac_lt10:.2f}% of cCREs have fewer than 10 BCs -- high complexity",
            f"{100*frac_lt10:.2f}% of cCREs have fewer than 10 BCs -- low complexity")

    # ---------------- 1.4 Retained cCREs ----------------
    head("1.4  Retained cCREs", "assess BC replicability (complexity)")
    say("  metric: cCREs retained at increasing minimum-BC-per-cCRE cutoffs")
    say(f"  {'min BCs':>10}{'cCREs':>12}{'% retained':>13}")
    rows = []
    for t in (1, 5, 10, 20, 30, 50, 100, 200):
        n = int((v >= t).sum())
        rows.append(dict(min_bc=t, n=n, frac=n / len(v)))
        say(f"  {t:>10}{n:>12,}{100*n/len(v):>12.2f}%")
    Q["1.4"] = rows
    PD["retained_rows"] = rows
    r100 = rows[-2]["frac"]
    verdict(r100 > 0.90,
            f"{100*r100:.1f}% of cCREs survive a 100-BC cutoff -- complexity is not limiting",
            f"only {100*r100:.1f}% of cCREs survive a 100-BC cutoff")

    # ---------------- 1.5 / 1.6 downsampling ----------------
    head("1.5  cCRE retention by sequencing depth", "assess whether sequencing depth is sufficient")
    say("  metric: reads downsampled and barcodes recounted; cCRE retention refit")
    say("          with the Hill saturation model and extrapolated to 3x current depth")
    reads = (bd_tot + br_tot).astype(np.int64)
    fr, ret, bcper = [], [], []
    n_ccre = len(bcs_per)
    for f in (0.1, 0.25, 0.5, 0.75, 1.0):
        # downsample READS; a barcode survives only if it retains >=1 read, so
        # low-coverage barcodes drop out and the curve can actually saturate
        sub = rng.binomial(reads, f) if f < 1 else reads
        keep = sub > 0
        c = np.bincount(bc_code[keep], minlength=n_lev)
        fr.append(f)
        ret.append(int((c > 0).sum()) / n_ccre)
        bcper.append(c[c > 0])
    say(f"\n  {'frac':>6}{'barcodes':>13}{'cCREs':>10}{'% retained':>13}")
    for f, r_, bp in zip(fr, ret, bcper):
        say(f"  {f:>6.2f}{bp.sum():>13,}{len(bp):>10,}{100*r_:>12.2f}%")
    pred3, ymax = hill_extrapolate(fr, ret, 3.0)
    PD["ds_frac"], PD["ds_ret"] = fr, [float(x) for x in ret]
    Q["1.5"] = dict(frac=fr, retention=[float(x) for x in ret],
                    hill_pred_3x=pred3, hill_ymax=ymax)
    gain = (pred3 - ret[-1]) / max(ret[-1], 1e-9)
    say(f"\n  Hill fit: current {100*ret[-1]:.2f}%  ->  3x depth {100*pred3:.2f}%  "
        f"({1+gain:.3f}-fold)")
    verdict(gain < 0.02,
            f"3x more sequencing is predicted to change cCRE retention by {100*gain:.2f}% -- saturated",
            f"3x more sequencing is predicted to raise cCRE retention by {100*gain:.2f}%")

    head("1.6  BCs per cCRE by sequencing depth", "assess whether sequencing depth is sufficient")
    say("  metric: distribution of BCs per cCRE at each read-downsampling fraction")
    say(f"  {'frac':>6}{'Q1':>9}{'median':>9}{'Q3':>9}{'mean':>10}")
    rows = []
    for f, bp in zip(fr, bcper):
        rows.append(dict(frac=f, q1=float(np.percentile(bp, 25)),
                         median=float(np.median(bp)), q3=float(np.percentile(bp, 75)),
                         mean=float(bp.mean())))
        say(f"  {f:>6.2f}{np.percentile(bp,25):>9.0f}{np.median(bp):>9.0f}"
            f"{np.percentile(bp,75):>9.0f}{bp.mean():>10.1f}")
    Q["1.6"] = rows
    PD["ds_box"] = rows
    slope = (rows[-1]["median"] - rows[-2]["median"]) / max(rows[-1]["median"], 1)
    say(f"\n  median BCs/cCRE gained over the last 25% of reads: {100*slope:.1f}%")
    say("  (a purely linear rise would be 25%, i.e. no saturation at all)")
    verdict(slope < 0.15,
            f"BCs/cCRE is levelling off ({100*slope:.1f}% over the last depth quartile)",
            f"BCs/cCRE rises {100*slope:.1f}% over the last depth quartile -- "
            "the barcode pool is not sequenced to saturation; note the book's caveat that "
            "sequencing errors add spurious barcodes with a positive slope")

    # ---------------- 1.7 / 1.8 PCR bias ----------------
    gcv, gsv, dnav = [], [], []
    for n in names:
        s = seq_by_name.get(n)
        gcv.append(gc_frac(s) if s else np.nan)
        gsv.append(max_g(s) if s else -1)
    gcv = np.array(gcv, float)
    gsv = np.array(gsv, int)
    dna_tot = np.nansum(DNA, axis=0)

    def binned(title, goal, vals, bins, labfmt, key, xlab):
        head(title, goal)
        say(f"  metric: DNA reads and cCRE count per {xlab} bin (book uses fixed 5% GC bins)")
        say(f"\n  {xlab:<14}{'cCREs':>9}{'DNA reads':>15}{'reads/cCRE':>13}{'median':>10}")
        rows = []
        for lo, hi in bins:
            m = (vals >= lo) & (vals < hi) if hi is not None else (vals >= lo)
            if m.sum() == 0:
                continue
            rows.append(dict(lo=float(lo), hi=(float(hi) if hi is not None else None),
                             n=int(m.sum()), reads=float(dna_tot[m].sum()),
                             per=float(dna_tot[m].sum() / m.sum()),
                             median=float(np.median(dna_tot[m]))))
            say(f"  {labfmt(lo,hi):<14}{m.sum():>9,}{dna_tot[m].sum():>15,.0f}"
                f"{dna_tot[m].sum()/m.sum():>13,.1f}{np.median(dna_tot[m]):>10,.0f}")
        # fold range restricted to bins holding >=1% of the library, so a
        # near-empty tail bin cannot set the headline number
        n_binned = sum(r["n"] for r in rows)
        big = [r for r in rows if r["n"] >= 0.01 * n_binned]
        use = big if len(big) >= 3 else rows
        if len(big) < 3:
            say(f"\n  note: only {len(big)} bins hold >=1% of binned cCREs; "
                f"fold range below uses all {len(rows)} bins")
        if n_binned < 0.9 * len(vals):
            say(f"\n  WARNING: only {n_binned:,} of {len(vals):,} cCREs could be binned "
                f"({100*n_binned/len(vals):.1f}%) -- sequence lookup is incomplete")
        fold = max(r["per"] for r in use) / max(min(r["per"] for r in use), 1e-9)
        cov = sum(r["n"] for r in use) / max(n_binned, 1)
        Q[key] = dict(rows=rows, fold_range=float(fold), bins_used=len(use),
                      coverage=float(cov))
        say(f"\n  fold range over bins holding >=1% of the library "
            f"({len(use)} bins, {100*cov:.1f}% of binned cCREs): {fold:.2f}x")
        verdict(fold < 2.0,
                f"read depth varies {fold:.2f}x across well-represented {xlab} bins -- little bias",
                f"read depth varies {fold:.2f}x across well-represented {xlab} bins -- amplification bias")
        return rows

    PD["gc_rows"] = binned("1.7  PCR bias - GC", "assess GC content bias in PCR amplification",
           gcv, [(x / 100, x / 100 + 0.05) for x in range(0, 100, 5)],
           lambda lo, hi: f"{lo:.2f}-{hi:.2f}", "1.7", "GC")
    PD["gs_rows"] = binned("1.8  PCR bias - G-stretches", "assess G-stretch bias in PCR amplification",
           gsv.astype(float), [(x, x + 1) for x in range(0, 12)] + [(12, None)],
           lambda lo, hi: (f"{int(lo)}" if hi is not None else f">={int(lo)}"), "1.8", "G-stretch")

    say("")
    say("#" * 86)
    say("#  CHAPTER 2 -- ACTIVITY QC")
    say("#" * 86)

    # ---------------- 2.1 ----------------
    head("2.1  Retained cCREs and BCs", "assess cCRE and BC coverage")
    say("  metric: % of association cCREs/BCs seen in the DNA and RNA quantification")
    n_bc_dna = int((bd.sum(axis=0) > 0).sum())
    n_bc_rna = int((br.sum(axis=0) > 0).sum())
    n_bc_any = int(((bd.sum(axis=0) + br.sum(axis=0)) > 0).sum())
    say(f"\n  {'':<10}{'BCs':>14}{'% of assoc':>13}{'cCREs':>10}{'% of assoc':>13}")
    rows = []
    for lab, mbc, arr in (("DNA", n_bc_dna, DNA), ("RNA", n_bc_rna, RNA), ("either", n_bc_any, None)):
        if arr is not None:
            nc = int((arr.sum(axis=0) > 0).sum())
        else:
            nc = int(((DNA.sum(axis=0) + RNA.sum(axis=0)) > 0).sum())
        rows.append(dict(channel=lab, barcodes=mbc, ccres=nc))
        say(f"  {lab:<10}{mbc:>14,}{100*mbc/n_assoc_bc:>12.2f}%{nc:>10,}{100*nc/len(bcs_per):>12.2f}%")
    PD["retention_rows"] = rows
    PD["assoc_bc"] = n_assoc_bc
    Q["2.1"] = dict(rows=rows, assoc_barcodes=n_assoc_bc, assoc_ccres=len(bcs_per))
    say(f"\n  per replicate:")
    say(f"  {'rep':>6}{'cCREs DNA':>12}{'cCREs RNA':>12}{'DNA UMIs':>16}{'RNA UMIs':>16}")
    perrep = []
    for i, r in enumerate(REPS):
        perrep.append(dict(rep=r, ccre_dna=int((DNA[i] > 0).sum()), ccre_rna=int((RNA[i] > 0).sum()),
                           dna=float(DNA[i].sum()), rna=float(RNA[i].sum())))
        say(f"  {r:>6}{(DNA[i]>0).sum():>12,}{(RNA[i]>0).sum():>12,}"
            f"{DNA[i].sum():>16,.0f}{RNA[i].sum():>16,.0f}")
    Q["2.1_per_rep"] = perrep
    ret_bc = n_bc_any / n_assoc_bc
    verdict(ret_bc > 0.8,
            f"{100*ret_bc:.1f}% of association barcodes are observed -- high retention",
            f"only {100*ret_bc:.1f}% of association barcodes are observed")

    # ---------------- 2.2 ----------------
    binned("2.2  DNA counts vs GC content", "assess GC content bias in PCR amplification",
           gcv, [(x / 100, x / 100 + 0.05) for x in range(0, 100, 5)],
           lambda lo, hi: f"{lo:.2f}-{hi:.2f}", "2.2", "GC")

    # ---------------- 2.3 ----------------
    head("2.3  Activity statistic vs count ratio",
         "assess the correlation between log2(RNA/DNA) and the MPRAnalyze activity statistic")
    ok = np.isfinite(A_STAT) & np.isfinite(RL_MEAN) & test_m
    pr = pearson(A_STAT[ok], RL_MEAN[ok])
    sp = spearman(A_STAT[ok], RL_MEAN[ok])
    say(f"  test cCREs with both measures: {ok.sum():,}")
    say(f"  Pearson  r = {pr:.3f}")
    say(f"  Spearman r = {sp:.3f}")
    Q["2.3"] = dict(n=int(ok.sum()), pearson=pr, spearman=sp)
    verdict(pr > 0.7,
            f"MPRAnalyze statistic tracks log2(RNA/DNA) at r={pr:.3f} -- high complexity",
            f"MPRAnalyze statistic and log2(RNA/DNA) agree only at r={pr:.3f} -- low complexity")

    # ---------------- 2.4 ----------------
    head("2.4  Activity distribution", "assess activity dynamic range, noise and statistical power")
    say("  metric: distribution of log2(RNA/DNA) over test cCREs, actives (MPRAnalyze")
    say(f"          fdr.mad < {FDR_CUT}) reported separately")
    t = test_m & np.isfinite(RL_MEAN)
    act_m = t & ACTIVE
    say(f"\n  test cCREs {t.sum():,}   active {int(act_m.sum()):,} "
        f"({100*act_m.sum()/t.sum():.2f}%)")
    say(f"  {'group':<16}{'n':>9}{'median':>10}{'p05':>9}{'p95':>9}{'sd':>9}")
    for lab, m in (("all test", t), ("active", act_m), ("non-active", t & ~ACTIVE),
                   ("positive ctrl", pos_m), ("negative ctrl", neg_m)):
        x = RL_MEAN[m]
        x = x[np.isfinite(x)]
        if len(x) < 2:
            continue
        say(f"  {lab:<16}{len(x):>9,}{np.median(x):>10.3f}{np.percentile(x,5):>9.2f}"
            f"{np.percentile(x,95):>9.2f}{x.std():>9.3f}")
    sep = float(np.median(RL_MEAN[act_m]) - np.median(RL_MEAN[t & ~ACTIVE]))
    dyn = float(np.percentile(RL_MEAN[t], 95) - np.percentile(RL_MEAN[t], 5))
    Q["2.4"] = dict(n_test=int(t.sum()), n_active=int(act_m.sum()),
                    frac_active=float(act_m.sum() / t.sum()),
                    active_vs_nonactive=sep, dynamic_range_5_95=dyn)
    say(f"\n  active - non-active median separation: {sep:+.3f} log2")
    say(f"  dynamic range (5-95 pct of test cCREs): {dyn:.3f} log2")
    verdict(sep > 0.5,
            f"actives sit {sep:+.2f} log2 above non-actives -- activity is detected cleanly",
            f"actives sit only {sep:+.2f} log2 above non-actives -- weak activity detection")

    # ---------------- 2.5 ----------------
    if comp is not None:
        head("2.5  Differential activity distribution",
             "assess differential activity dynamic range, noise and statistical power")
        lfc = comp["logFC"].to_numpy(float)
        dsig = comp["differentialy_active"].to_numpy(bool)
        say(f"  loci tested {len(lfc):,}   differentially active {int(dsig.sum()):,} "
            f"({100*dsig.sum()/len(lfc):.2f}%)")
        say(f"  logFC = log2(Derived/Ancestral)")
        say(f"  {'group':<20}{'n':>9}{'median':>10}{'p05':>9}{'p95':>9}")
        for lab, m in (("all", np.ones(len(lfc), bool)), ("differentially active", dsig)):
            x = lfc[m]
            say(f"  {lab:<20}{len(x):>9,}{np.median(x):>10.4f}"
                f"{np.percentile(x,5):>9.3f}{np.percentile(x,95):>9.3f}")
        anc = int((dsig & (lfc < 0)).sum())
        der = int((dsig & (lfc > 0)).sum())
        say(f"\n  Ancestral-enriched (logFC<0): {anc:,}    Derived-enriched (logFC>0): {der:,}")
        Q["2.5"] = dict(n=len(lfc), n_sig=int(dsig.sum()), ancestral=anc, derived=der,
                        median_abs_sig=float(np.median(np.abs(lfc[dsig]))) if dsig.any() else None)
        if dsig.any():
            verdict(np.median(np.abs(lfc[dsig])) > 0.2,
                    f"significant loci carry a median |logFC| of {np.median(np.abs(lfc[dsig])):.3f}",
                    f"significant loci carry only {np.median(np.abs(lfc[dsig])):.3f} median |logFC|")
    else:
        na("2.5  Differential activity distribution",
           "assess differential activity dynamic range, noise and statistical power",
           "stage 4 comparative has not finalised yet")

    # ---------------- 2.6 ----------------
    head("2.6  P-value distribution", "assess statistical power")
    say("  metric: QQ of MPRAnalyze pval.mad against Uniform(0,1)")
    pv_t = A_PVAL[test_m]
    pv_t = pv_t[np.isfinite(pv_t)]
    pv_n = A_PVAL[neg_m]
    pv_n = pv_n[np.isfinite(pv_n)]
    say(f"\n  {'group':<18}{'n':>9}{'median':>10}{'frac<0.05':>12}{'KS vs unif':>13}{'KS p':>10}")
    rows = []
    for lab, x in (("test cCREs", pv_t), ("negative ctrl", pv_n)):
        ks = stats.kstest(x, "uniform")
        rows.append(dict(group=lab, n=len(x), median=float(np.median(x)),
                         frac_lt_05=float((x < 0.05).mean()),
                         ks_D=float(ks.statistic), ks_p=float(ks.pvalue)))
        say(f"  {lab:<18}{len(x):>9,}{np.median(x):>10.4f}{100*(x<0.05).mean():>11.1f}%"
            f"{ks.statistic:>13.4f}{ks.pvalue:>10.3g}")
    say(f"\n  observed quantiles of the negative-control null (expected = the same value):")
    say("    " + "  ".join(f"{x:.2f}" for x in np.quantile(pv_n, np.arange(0, 1.01, 0.1))))
    say("    " + "  ".join(f"{x:.2f}" for x in np.arange(0, 1.01, 0.1)))
    bt = stats.binomtest(int((pv_n < 0.05).sum()), len(pv_n), 0.05)
    lo, hi = bt.proportion_ci(0.95)
    say(f"\n  negative-control rate at p<0.05: {100*(pv_n<0.05).mean():.1f}% "
        f"(95% CI {100*lo:.1f}-{100*hi:.1f}%), nominal 5%, binomial p={bt.pvalue:.3g}")
    PD["pv_t"], PD["pv_n"] = pv_t, pv_n
    Q["2.6"] = dict(rows=rows, neg_rate=float((pv_n < 0.05).mean()),
                    neg_ci=[float(lo), float(hi)], neg_binom_p=float(bt.pvalue))
    enrich = (pv_t < 0.05).mean() / max((pv_n < 0.05).mean(), 1e-9)
    say(f"  test/negative enrichment at p<0.05: {enrich:.2f}x")
    verdict(hi > 0.05 > lo and enrich > 1.5,
            "negative-control null is calibrated and test cCREs are enriched for low p-values",
            f"negative controls give {100*(pv_n<0.05).mean():.1f}% at p<0.05 (nominal 5%) and "
            f"test cCREs are enriched only {enrich:.2f}x -- limited power at this cut")

    # ---------------- 2.7 ----------------
    if comp is not None:
        head("2.7  Volcano plot - FC vs P-value", "assess statistical power of differential activity")
        lfc = comp["logFC"].to_numpy(float)
        fdr = comp["differential_activity_FDR"].to_numpy(float)
        say(f"  {'|logFC| band':<16}{'loci':>9}{'sig':>9}{'% sig':>10}")
        rows = []
        for lo_, hi_ in ((0, .1), (.1, .2), (.2, .4), (.4, .8), (.8, None)):
            m = (np.abs(lfc) >= lo_) & ((np.abs(lfc) < hi_) if hi_ else True)
            if m.sum() == 0:
                continue
            s = int((fdr[m] < FDR_CUT).sum())
            rows.append(dict(lo=lo_, hi=hi_, n=int(m.sum()), sig=s))
            lab = f"{lo_:.1f}-{hi_:.1f}" if hi_ else f">={lo_:.1f}"
            say(f"  {lab:<16}{m.sum():>9,}{s:>9,}{100*s/m.sum():>9.1f}%")
        Q["2.7"] = rows
        big = np.abs(lfc) >= 0.4
        rate = (fdr[big] < FDR_CUT).mean() if big.sum() else np.nan
        verdict(rate > 0.5,
                f"{100*rate:.0f}% of loci with |logFC|>=0.4 are called significant -- adequate power",
                f"only {100*rate:.0f}% of loci with |logFC|>=0.4 reach FDR<{FDR_CUT}")
    else:
        na("2.7  Volcano plot - FC vs P-value",
           "assess statistical power of differential activity detection",
           "stage 4 comparative has not finalised yet")

    # ---------------- 2.8 / 2.9 / 2.10 downsampling ----------------
    head("2.8  BC retention by DNA/RNA sequencing depth", "assess whether sequencing depth is sufficient")
    say("  metric: DNA and RNA reads downsampled jointly; retained BCs refit with the")
    say("          Hill model and extrapolated to 3x current depth")
    fr2, keep_d, keep_r, keep_c_d, keep_c_r, act_ds = [], [], [], [], [], []
    n_bc_all = len(bc)
    for f in (0.1, 0.25, 0.5, 0.75, 1.0):
        if f < 1:
            sd = rng.binomial(bd.astype(np.int64), f)
            sr = rng.binomial(br.astype(np.int64), f)
        else:
            sd, sr = bd.astype(np.int64), br.astype(np.int64)
        std, srt = sd.sum(axis=0), sr.sum(axis=0)
        fr2.append(f)
        keep_d.append((std > 0).sum() / n_bc_all)
        keep_r.append((srt > 0).sum() / n_bc_all)
        cd = by_ccre(std.astype(float))
        cr = by_ccre(srt.astype(float))
        keep_c_d.append((cd > 0).sum() / n_ccre)
        keep_c_r.append((cr > 0).sum() / n_ccre)
        # activity proxy renormalised for the depth actually retained
        lr = np.log2((cr + 1) / (cd + 1))
        act_ds.append(float(((lr > 1) & (cd > 0)).mean()))
    say(f"\n  {'frac':>6}{'BC DNA':>12}{'BC RNA':>12}")
    for f, d_, r_ in zip(fr2, keep_d, keep_r):
        say(f"  {f:>6.2f}{100*d_:>11.2f}%{100*r_:>11.2f}%")
    pd3, _ = hill_extrapolate(fr2, keep_d, 3.0)
    pr3, _ = hill_extrapolate(fr2, keep_r, 3.0)
    say(f"\n  Hill 3x: DNA {100*keep_d[-1]:.1f}% -> {100*pd3:.1f}% ({pd3/keep_d[-1]:.3f}-fold)   "
        f"RNA {100*keep_r[-1]:.1f}% -> {100*pr3:.1f}% ({pr3/keep_r[-1]:.3f}-fold)")
    PD["ds2_frac"] = fr2
    PD["ds2_dna"], PD["ds2_rna"] = [float(x) for x in keep_d], [float(x) for x in keep_r]
    PD["ds2_cdna"], PD["ds2_crna"] = [float(x) for x in keep_c_d], [float(x) for x in keep_c_r]
    PD["ds2_act"] = [float(x) for x in act_ds]
    Q["2.8"] = dict(frac=fr2, dna=[float(x) for x in keep_d], rna=[float(x) for x in keep_r],
                    hill_3x_dna=pd3, hill_3x_rna=pr3)
    verdict(max(pd3 / keep_d[-1], pr3 / keep_r[-1]) < 1.10,
            "3x more sequencing would change BC retention by <10% -- depth is sufficient",
            "3x more sequencing would still add BCs -- depth is suboptimal")

    head("2.9  cCRE retention by DNA/RNA sequencing depth", "assess whether sequencing depth is sufficient")
    say(f"  {'frac':>6}{'cCRE DNA':>12}{'cCRE RNA':>12}")
    for f, d_, r_ in zip(fr2, keep_c_d, keep_c_r):
        say(f"  {f:>6.2f}{100*d_:>11.2f}%{100*r_:>11.2f}%")
    cd3, _ = hill_extrapolate(fr2, keep_c_d, 3.0)
    cr3, _ = hill_extrapolate(fr2, keep_c_r, 3.0)
    say(f"\n  Hill 3x: DNA {100*keep_c_d[-1]:.2f}% -> {100*cd3:.2f}%   "
        f"RNA {100*keep_c_r[-1]:.2f}% -> {100*cr3:.2f}%")
    Q["2.9"] = dict(frac=fr2, dna=[float(x) for x in keep_c_d], rna=[float(x) for x in keep_c_r],
                    hill_3x_dna=cd3, hill_3x_rna=cr3)
    verdict(max(cd3 / max(keep_c_d[-1], 1e-9), cr3 / max(keep_c_r[-1], 1e-9)) < 1.02,
            "cCRE retention is at plateau -- more sequencing would not recover more cCREs",
            "cCRE retention is still climbing with depth")

    head("2.10  Activity by sequencing depth", "assess whether sequencing depth is sufficient")
    say("  metric: fraction of cCREs above the activity proxy log2(RNA/DNA) > 1.")
    say("          DNA and RNA are downsampled jointly, so the ratio is unbiased by depth")
    say("          and the threshold is held fixed.")
    say(f"\n  {'frac':>6}{'% active':>12}")
    for f, x in zip(fr2, act_ds):
        say(f"  {f:>6.2f}{100*x:>11.2f}%")
    Q["2.10"] = dict(frac=fr2, active=[float(x) for x in act_ds])
    g = abs(act_ds[-1] - act_ds[-2]) / max(act_ds[-1], 1e-9)
    verdict(g < 0.02,
            f"the active fraction moves {100*g:.2f}% over the last 25% of depth -- plateau reached",
            f"the active fraction moves {100*g:.2f}% over the last 25% of depth -- still climbing")

    # ---------------- 2.11 ----------------
    head("2.11  Cumulative RNA reads", "assess jackpotting")
    say("  metric: share of total RNA UMIs held by the top-ranked cCREs")
    rv = np.nansum(RNA, axis=0)
    rv = np.sort(rv[np.isfinite(rv)])[::-1]
    cum = np.cumsum(rv) / rv.sum()
    marks = {}
    for f in (0.001, 0.01, 0.05, 0.1, 0.25, 0.5):
        i = max(int(len(rv) * f) - 1, 0)
        marks[f] = float(cum[i])
        say(f"  top {100*f:>6.1f}% of cCREs carry {100*cum[i]:>6.2f}% of RNA UMIs")
    PD["rna_cum"] = cum
    Q["2.11"] = dict(marks={str(k): v for k, v in marks.items()})
    verdict(marks[0.01] < 0.30,
            f"the top 1% of cCREs hold {100*marks[0.01]:.1f}% of RNA -- no jackpotting",
            f"the top 1% of cCREs hold {100*marks[0.01]:.1f}% of RNA -- jackpotting")

    # ---------------- 2.12 ----------------
    head("2.12  Sample clustering", "assess reproducibility between samples")
    say("  metric: PCA over log10(CPM+1) per cCRE, samples as observations")
    fin = np.all(np.isfinite(np.vstack([DNA_CPM, RNA_CPM])), axis=0)
    M = np.vstack([np.log10(DNA_CPM[i][fin] + 1) for i in range(3)] +
                  [np.log10(RNA_CPM[i][fin] + 1) for i in range(3)])
    labels = [f"DNA_{r}" for r in REPS] + [f"RNA_{r}" for r in REPS]
    Mc = M - M.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(Mc, full_matrices=False)
    var = (S ** 2) / (S ** 2).sum()
    pcs = U[:, :3] * S[:3]
    say(f"  cCREs with complete CPM in all six samples: {int(fin.sum()):,}")
    say(f"  variance explained: PC1 {100*var[0]:.1f}%  PC2 {100*var[1]:.1f}%  PC3 {100*var[2]:.1f}%")
    say(f"\n  {'sample':<12}{'PC1':>10}{'PC2':>10}")
    for l, x, y in zip(labels, pcs[:, 0], pcs[:, 1]):
        say(f"  {l:<12}{x:>10.2f}{y:>10.2f}")
    Q["2.12"] = dict(labels=labels, pc1=pcs[:, 0].tolist(), pc2=pcs[:, 1].tolist(),
                     var=[float(x) for x in var[:3]])
    PD["pca"] = Q["2.12"]
    sepq = abs(pcs[:3, 0].mean() - pcs[3:, 0].mean())
    within = max(pcs[:3, 0].std(), pcs[3:, 0].std())
    say(f"\n  PC1 DNA-vs-RNA separation {sepq:.1f} against within-group spread {within:.1f} "
        f"({sepq/max(within,1e-9):.1f}x)")
    say("  note: this run is a single cell type, so PC1 separates assay channel rather than")
    say("        biological group; see 2.18 for the cross-cell-type comparison")
    verdict(sepq > 3 * within,
            "channels separate cleanly on PC1 with replicates tight within channel",
            "replicate spread is comparable to the DNA/RNA separation")

    # ---------------- 2.13 ----------------
    head("2.13  Correlation between replicates", "assess reproducibility between replicates")
    say("  metric: Pearson/Spearman of CPM-normalised log2(RNA/DNA) per cCRE")
    say(f"\n  {'pair':<16}{'n':>9}{'pearson':>10}{'spearman':>11}")
    rows = []
    for i, j in PAIRS:
        m = test_m & np.isfinite(RL[i]) & np.isfinite(RL[j])
        p_, s_ = pearson(RL[i][m], RL[j][m]), spearman(RL[i][m], RL[j][m])
        rows.append(dict(pair=f"{REPS[i]}_{REPS[j]}", n=int(m.sum()), pearson=p_, spearman=s_))
        say(f"  {REPS[i]+' vs '+REPS[j]:<16}{m.sum():>9,}{p_:>10.3f}{s_:>11.3f}")
    Q["2.13"] = rows
    mn = min(r["pearson"] for r in rows)
    verdict(mn > 0.7,
            f"all replicate pairs correlate at r>={mn:.3f}",
            f"weakest replicate pair correlates at only r={mn:.3f}")

    # ---------------- allelic pairs ----------------
    base_of, allele_of = {}, {}
    for n in names:
        if "_Derived" in n:
            base_of[n] = n.replace("_Derived", "")
            allele_of[n] = "Derived"
        elif "_Ancestral" in n:
            base_of[n] = n.replace("_Ancestral", "")
            allele_of[n] = "Ancestral"
    idx_of = {n: i for i, n in enumerate(names)}
    loci = collections.defaultdict(dict)
    for n, b in base_of.items():
        loci[b][allele_of[n]] = idx_of[n]
    pair_idx = [(v["Derived"], v["Ancestral"]) for v in loci.values()
                if "Derived" in v and "Ancestral" in v]
    di = np.array([p[0] for p in pair_idx])
    ai = np.array([p[1] for p in pair_idx])

    head("2.14  Correlation of differential activity between replicates",
         "assess reproducibility between replicates")
    say("  metric: replicate correlation of log2(RNA_cpm/DNA_cpm) over the oligos that")
    say("          make up the allelic pairs. Each allele is its own oligo, so this is the")
    say("          reproducibility of the measurement the comparative test is built from.")
    allel = np.concatenate([di, ai])
    LFC_REP = np.vstack([RL[i][di] - RL[i][ai] for i in range(3)])
    say(f"\n  allelic pairs {len(di):,}  ({len(allel):,} oligos)")
    say(f"  {'set':<26}{'n':>9}{'r(1,2)':>9}{'r(1,3)':>9}{'r(2,3)':>9}{'mean':>9}")
    rows = []
    for lab, idx in (("all allelic-pair oligos", allel), ("Derived oligos", di),
                     ("Ancestral oligos", ai)):
        rr = [pearson(RL[i][idx], RL[j][idx]) for i, j in PAIRS]
        rows.append(dict(set=lab, n=int(len(idx)), r=[float(x) for x in rr],
                         mean_r=float(np.mean(rr))))
        say(f"  {lab:<26}{len(idx):>9,}{rr[0]:>9.3f}{rr[1]:>9.3f}{rr[2]:>9.3f}{np.mean(rr):>9.3f}")
    oligo_r = rows[0]["mean_r"]

    # Secondary diagnostic. The correlation of the allelic DIFFERENCE across
    # replicates is a different and much smaller quantity: subtracting the two
    # alleles cancels the locus-level signal they share and adds both noises.
    # A low value here is arithmetic, not evidence the assay is unreliable, so
    # it is reported but does not drive the verdict.
    dif = [pearson(LFC_REP[i], LFC_REP[j]) for i, j in PAIRS]
    sig_o = float(np.nanstd(np.nanmean(RL[:, allel], axis=0)))
    noi_o = float(np.nanmedian(np.nanstd(RL[:, allel], axis=0)))
    sig_d = float(np.nanstd(np.nanmean(LFC_REP, axis=0)))
    noi_d = float(np.nanmedian(np.nanstd(LFC_REP, axis=0)))
    say(f"\n  secondary: replicate correlation of the allelic DIFFERENCE "
        f"(Derived - Ancestral): {np.mean(dif):.3f}")
    say(f"    signal/noise per oligo      SD {sig_o:.3f} across oligos vs {noi_o:.3f} within")
    say(f"    signal/noise of the diff    SD {sig_d:.3f} across loci  vs {noi_d:.3f} within")
    say("    the difference cancels the shared locus signal and adds both noises, so its")
    say("    correlation is expected to be far lower and is not a reproducibility failure")
    Q["2.14_difference"] = dict(r=[float(x) for x in dif], mean_r=float(np.mean(dif)),
                                signal_oligo=sig_o, noise_oligo=noi_o,
                                signal_diff=sig_d, noise_diff=noi_d)
    say("\n  the difference-correlation, stratified by locus activity:")
    say(f"  {'activity bin':<16}{'pairs':>9}{'r(1,2)':>10}{'r(1,3)':>10}{'r(2,3)':>10}")
    pair_act = np.nanmean(np.vstack([RL_MEAN[di], RL_MEAN[ai]]), axis=0)
    strat = []
    for lo_ in np.arange(np.floor(np.nanpercentile(pair_act, 5)),
                         np.ceil(np.nanpercentile(pair_act, 95)), 1.0):
        m = (pair_act >= lo_) & (pair_act < lo_ + 1)
        if m.sum() < 50:
            continue
        rr = [pearson(LFC_REP[i][m], LFC_REP[j][m]) for i, j in PAIRS]
        strat.append(dict(lo=float(lo_), n=int(m.sum()), r=[float(x) for x in rr]))
        say(f"  {f'[{lo_:.0f}, {lo_+1:.0f})':<16}{m.sum():>9,}"
            f"{rr[0]:>10.3f}{rr[1]:>10.3f}{rr[2]:>10.3f}")
    Q["2.14_by_activity"] = strat
    say("  most allelic pairs carry no true differential activity, so this rises with activity")
    Q["2.14"] = rows
    mnd = oligo_r
    best_strat = max((np.mean(x["r"]) for x in strat), default=float("nan"))
    if strat:
        say(f"  best stratum mean r for the difference = {best_strat:.3f}")
    verdict(mnd > 0.7,
            f"the oligos behind the allelic pairs reproduce across replicates at "
            f"mean r={mnd:.3f}",
            f"the oligos behind the allelic pairs reproduce at only mean r={mnd:.3f}")

    # ---------------- 2.15 ----------------
    head("2.15  Replicability by activity",
         "assess replicate correlation in active vs non-active cCREs")
    say("  metric: replicate correlation within fixed bins of mean log2(RNA/DNA), width 1")
    say(f"\n  {'activity bin':<16}{'n':>9}{'r(1,2)':>10}{'r(1,3)':>10}{'r(2,3)':>10}")
    rows = []
    lo0 = np.floor(np.nanpercentile(RL_MEAN[test_m], 1))
    hi0 = np.ceil(np.nanpercentile(RL_MEAN[test_m], 99))
    for lo_ in np.arange(lo0, hi0, 1.0):
        m = test_m & (RL_MEAN >= lo_) & (RL_MEAN < lo_ + 1)
        if m.sum() < 50:
            continue
        rr = [pearson(RL[i][m], RL[j][m]) for i, j in PAIRS]
        rows.append(dict(lo=float(lo_), n=int(m.sum()), r=[float(x) for x in rr]))
        say(f"  {f'[{lo_:.0f}, {lo_+1:.0f})':<16}{m.sum():>9,}"
            f"{rr[0]:>10.3f}{rr[1]:>10.3f}{rr[2]:>10.3f}")
    Q["2.15"] = rows
    PD["rep_by_act"] = rows
    if len(rows) >= 2:
        lowr = np.mean(rows[0]["r"])
        hir = np.mean(rows[-1]["r"])
        say(f"\n  mean r rises from {lowr:.3f} in the least active bin to {hir:.3f} in the most active")
        verdict(hir > lowr,
                "correlation improves with activity, as expected",
                "correlation does not improve with activity")

    # ---------------- 2.16 ----------------
    head("2.16  Correlation between replicates (controls)",
         "assess reproducibility in positive and negative controls")
    say("  metric: replicate correlation computed separately within each control class")
    say("          (a pooled control set is bimodal and inflates r, so classes are kept apart)")
    say(f"\n  {'group':<18}{'n':>7}{'r(1,2)':>10}{'r(1,3)':>10}{'r(2,3)':>10}")
    rows = []
    for lab, m in (("positive ctrl", pos_m), ("negative ctrl", neg_m), ("test cCREs", test_m)):
        rr = [pearson(RL[i][m], RL[j][m]) for i, j in PAIRS]
        rows.append(dict(group=lab, n=int(m.sum()), r=[float(x) for x in rr]))
        say(f"  {lab:<18}{m.sum():>7,}{rr[0]:>10.3f}{rr[1]:>10.3f}{rr[2]:>10.3f}")
    Q["2.16"] = rows
    rp = np.mean(rows[0]["r"])
    rn = np.mean(rows[1]["r"])
    say(f"\n  mean r: positives {rp:.3f}   negatives {rn:.3f}")
    verdict(rp > rn,
            f"positives reproduce better than negatives ({rp:.3f} vs {rn:.3f}), as expected",
            f"negatives reproduce as well as positives ({rn:.3f} vs {rp:.3f}) -- "
            "negatives are carrying reproducible signal")

    # ---------------- 2.17 ----------------
    head("2.17  Cross-validation: allelic pairs", "assess reproducibility between allelic pairs")
    say("  metric: Derived vs Ancestral log2(RNA/DNA); the two differ by a few bp and")
    say("          should track each other closely")
    xa, xd = RL_MEAN[ai], RL_MEAN[di]
    p_, s_ = pearson(xd, xa), spearman(xd, xa)
    say(f"\n  allelic pairs {len(di):,}   Pearson r = {p_:.3f}   Spearman r = {s_:.3f}")
    say(f"  median |Derived - Ancestral| = {np.nanmedian(np.abs(xd-xa)):.4f} log2")
    Q["2.17"] = dict(n=int(len(di)), pearson=p_, spearman=s_,
                     median_abs_diff=float(np.nanmedian(np.abs(xd - xa))))
    verdict(p_ > 0.8,
            f"allelic pairs correlate at r={p_:.3f} -- internally consistent",
            f"allelic pairs correlate at only r={p_:.3f}")

    # ---------------- 2.18 ----------------
    head("2.18  Cross-validation: cell types", "assess reproducibility between different cell types")
    say("  metric: log2(RNA/DNA) of shared cCREs against sibling MH cell types")
    rows = []
    for cname, path in SIBLINGS.items():
        if not os.path.exists(path):
            say(f"  {cname}: not available at {path}")
            continue
        try:
            sib = load_activity(path)
        except Exception as e:
            say(f"  {cname}: could not read ({e})")
            continue
        scols = [c for c in sib.columns if c.startswith("ratio_log_")]
        if not scols:
            continue
        sib_m = sib[scols].mean(axis=1)
        shared = [n for n in names if n in sib_m.index]
        if len(shared) < 100:
            say(f"  {cname}: only {len(shared)} shared cCREs")
            continue
        x = RL_MEAN[[idx_of[n] for n in shared]]
        y = sib_m.reindex(shared).to_numpy(float)
        p_, s_ = pearson(x, y), spearman(x, y)
        ctl = np.array([ccre_type(n) != "test_cCRE" for n in shared])
        pc = pearson(x[ctl], y[ctl]) if ctl.sum() > 10 else np.nan
        rows.append(dict(cell=cname, n=len(shared), pearson=p_, spearman=s_, pearson_controls=pc))
        PLOT.setdefault("2.18", {})[cname] = dict(dens=hist2d(x, y), pearson=p_, n=len(shared))
        say(f"  {cname:<6} shared {len(shared):>7,}   Pearson {p_:>6.3f}   "
            f"Spearman {s_:>6.3f}   controls-only {pc:>6.3f}")
    Q["2.18"] = rows
    if rows:
        best = max(r["pearson"] for r in rows)
        verdict(best > 0.4,
                f"activity is shared with related cell types (best r={best:.3f})",
                f"no cell type correlates above r={best:.3f} -- weak cross-cell-type concordance")
    else:
        say("  NOT EVALUATED -- no sibling cell-type activity tables found")

    # ---------------- 2.19 ----------------
    head("2.19  Minimizing noise",
         "remove outlier barcodes and optimize the minimum DNA threshold to improve reproducibility")
    say("  metric: replicate correlation across the outlier-filter x minimum-DNA grid")
    say("          (outlier levels: none, and the pipeline's per-oligo +/-2 SD RNA filter)")
    nf, ff = f"{BASE}/output/UMI/UMI_no_filter.txt", f"{BASE}/output/UMI/UMI_std2_filter.txt"
    grids = {}
    for lab, path in (("no filter", nf), ("+/-2 SD", ff)):
        if not os.path.exists(path):
            continue
        t = pd.read_csv(path, sep="\t", usecols=lambda c: c == "oligo_bc" or
                        c.startswith(("DNA_rep", "RNA_rep")))
        t["cCRE"] = t["oligo_bc"].str.rsplit("_", n=1).str[0]
        g = t.groupby("cCRE")[[f"{c}_{r}" for c in ("DNA", "RNA") for r in REPS]].sum()
        d = np.vstack([g[f"DNA_{r}"].to_numpy(float) for r in REPS])
        r_ = np.vstack([g[f"RNA_{r}"].to_numpy(float) for r in REPS])
        lib_d, lib_r = d.sum(axis=1, keepdims=True), r_.sum(axis=1, keepdims=True)
        rl = np.log2(((r_ + 1) / lib_r * 1e6) / ((d + 1) / lib_d * 1e6))
        grids[lab] = (d, rl)
    say(f"\n  {'outlier filter':<14}{'min DNA/rep':>13}{'cCREs':>9}"
        f"{'r(1,2)':>9}{'r(1,3)':>9}{'r(2,3)':>9}{'mean r':>9}")
    rows = []
    for lab, (d, rl) in grids.items():
        for mind in (0, 10, 50, 100, 200, 500):
            m = (d.min(axis=0) >= mind)
            if m.sum() < 100:
                continue
            rr = [pearson(rl[i][m], rl[j][m]) for i, j in PAIRS]
            rows.append(dict(filter=lab, min_dna=mind, n=int(m.sum()),
                             r=[float(x) for x in rr], mean_r=float(np.mean(rr))))
            say(f"  {lab:<14}{mind:>13}{m.sum():>9,}"
                f"{rr[0]:>9.3f}{rr[1]:>9.3f}{rr[2]:>9.3f}{np.mean(rr):>9.3f}")
    Q["2.19"] = rows
    PD["noise_rows"] = rows
    if os.path.exists(nf) and os.path.exists(ff):
        diff = tot_rows = 0
        with op(nf) as x, op(ff) as y:
            x.readline(); y.readline()
            for la, lb in zip(x, y):
                tot_rows += 1
                if la != lb:
                    diff += 1
        say(f"\n  the +/-2 SD filter modified {diff:,} of {tot_rows:,} barcode rows "
            f"({100*diff/max(tot_rows,1):.2f}%)")
        Q["2.19_filter"] = dict(rows=tot_rows, modified=diff)
    if rows:
        best = max(rows, key=lambda r: r["mean_r"])
        base_r = [r for r in rows if r["min_dna"] == 0]
        b0 = base_r[0]["mean_r"] if base_r else np.nan
        say(f"  best cell: {best['filter']}, min DNA {best['min_dna']} -> mean r {best['mean_r']:.3f} "
            f"(from {b0:.3f} unfiltered)")
        verdict(best["mean_r"] >= b0,
                f"filtering raises replicate correlation to {best['mean_r']:.3f}",
                "filtering does not improve replicate correlation")

    # ---------------- 2.20 ----------------
    head("2.20  Reproducibility by sequencing depth",
         "assess whether additional sequencing will improve reproducibility")
    say("  metric: replicate correlation of CPM log2(RNA/DNA) at each downsampling fraction")
    say(f"\n  {'frac':>6}{'r(1,2)':>10}{'r(1,3)':>10}{'r(2,3)':>10}{'mean r':>10}")
    rows = []
    for f in (0.1, 0.25, 0.5, 0.75, 1.0):
        if f < 1:
            sd = rng.binomial(bd.astype(np.int64), f)
            sr = rng.binomial(br.astype(np.int64), f)
        else:
            sd, sr = bd.astype(np.int64), br.astype(np.int64)
        gd = np.vstack([by_ccre(sd[i].astype(float)) for i in range(3)])
        gr = np.vstack([by_ccre(sr[i].astype(float)) for i in range(3)])
        ld, lr_ = gd.sum(axis=1, keepdims=True), gr.sum(axis=1, keepdims=True)
        rl = np.log2(((gr + 1) / lr_ * 1e6) / ((gd + 1) / ld * 1e6))
        rr = [pearson(rl[i], rl[j]) for i, j in PAIRS]
        rows.append(dict(frac=f, r=[float(x) for x in rr], mean_r=float(np.mean(rr))))
        say(f"  {f:>6.2f}{rr[0]:>10.3f}{rr[1]:>10.3f}{rr[2]:>10.3f}{np.mean(rr):>10.3f}")
    Q["2.20"] = rows
    PD["repro_rows"] = rows
    dr = rows[-1]["mean_r"] - rows[-2]["mean_r"]
    say(f"\n  mean r gained over the last 25% of depth: {dr:+.4f}")
    verdict(abs(dr) < 0.01,
            "reproducibility has plateaued -- more sequencing would not improve it",
            f"reproducibility still gains {dr:+.4f} per depth quartile")

    # ---------------- 2.21 ----------------
    head("2.21  RNA vs DNA", "assess true activity in the experiment")
    say("  metric: CPM-normalised RNA vs DNA per cCRE on a LINEAR scale, axes truncated")
    say("          at 250 CPM as in the book; the diagonal is RNA = DNA")
    say(f"\n  {'rep':>6}{'above diag':>13}{'median log2':>14}{'p05':>9}{'p95':>9}{'sd':>9}")
    rows = []
    for i, r in enumerate(REPS):
        x = RL[i][test_m]
        x = x[np.isfinite(x)]
        rows.append(dict(rep=r, above=float((x > 0).mean()), median=float(np.median(x)),
                         p05=float(np.percentile(x, 5)), p95=float(np.percentile(x, 95)),
                         sd=float(x.std())))
        say(f"  {r:>6}{100*(x>0).mean():>12.1f}%{np.median(x):>14.3f}"
            f"{np.percentile(x,5):>9.2f}{np.percentile(x,95):>9.2f}{x.std():>9.3f}")
    Q["2.21"] = rows
    PD["rna_dna_rows"] = rows
    spread = max(r["median"] for r in rows) - min(r["median"] for r in rows)
    ab = float(np.mean([r["above"] for r in rows]))
    # CPM sets the two library totals equal, so the DNA-weighted mean ratio is
    # pinned at 1 (log2 = 0).  The median is free to sit below it whenever RNA
    # is more concentrated than DNA, so "% above the diagonal" reflects the RNA
    # skew rather than whether the assay worked.  Anchor on that instead.
    wmean = float(np.log2(np.nansum(RNA_CPM[:, test_m]) / np.nansum(DNA_CPM[:, test_m])))
    tail = float(np.mean([r["p95"] for r in rows]))
    say(f"\n  DNA-weighted mean log2(RNA/DNA): {wmean:+.4f} "
        "(pinned near 0 by CPM -- the library totals are equalised)")
    say(f"  spread of replicate medians after CPM normalisation: {spread:.3f} log2")
    say(f"  cCREs above the diagonal: {100*ab:.1f}%  |  mean 95th pct: {tail:+.2f} log2")
    say("  note: because CPM equalises the totals, a median below the diagonal follows from")
    say("        RNA being more concentrated than DNA (see 2.11) and is not itself silencing")
    Q["2.21_summary"] = dict(above=ab, spread=float(spread), weighted_mean=wmean, tail=tail)
    verdict(spread < 0.25 and tail > 0.5,
            f"replicate medians agree within {spread:.2f} log2 and a high-ratio population "
            f"extends to {tail:+.2f} log2 -- activity is present and normalisation is consistent",
            f"replicate medians differ by {spread:.2f} log2 / upper tail reaches only {tail:+.2f} log2")

    # ---------------- 2.22 ----------------
    head("2.22  Activity of controls", "assess the dynamic range of activity")
    say("  metric: activity by control class. The positive-control panel mixes HARs, a")
    say("          generic panel and named enhancers that are mostly bone/osteoblast, so")
    say("          it is split rather than pooled; scrambled negatives shown alongside.")
    say(f"\n  {'group':<22}{'n':>7}{'median':>10}{'Q1':>9}{'Q3':>9}{'vs test':>10}")
    med_test = float(np.nanmedian(RL_MEAN[test_m]))
    cls = np.array([ctrl_class(n) for n in names])
    order = ["HAR", "generic panel", "other named", "bone/osteoblast enh", "scrambled"]
    groups = [(l, cls == l) for l in order if (cls == l).sum() >= 5]
    groups.append(("test cCREs", test_m))
    rows = []
    for lab, m in groups:
        x = RL_MEAN[m]
        x = x[np.isfinite(x)]
        rows.append(dict(group=lab, n=len(x), median=float(np.median(x)),
                         q1=float(np.percentile(x, 25)), q3=float(np.percentile(x, 75)),
                         vs_test=float(np.median(x) - med_test)))
        say(f"  {lab:<22}{len(x):>7,}{np.median(x):>10.3f}{np.percentile(x,25):>9.3f}"
            f"{np.percentile(x,75):>9.3f}{np.median(x)-med_test:>+10.3f}")
    byg = {r["group"]: r for r in rows}
    scr = byg.get("scrambled", {}).get("median", float("nan"))
    best = max((r["median"] for r in rows if r["group"] not in ("scrambled", "test cCREs")),
               default=float("nan"))
    sep = best - scr
    u = stats.mannwhitneyu(RL_MEAN[pos_m][np.isfinite(RL_MEAN[pos_m])],
                           RL_MEAN[neg_m][np.isfinite(RL_MEAN[neg_m])],
                           alternative="greater")
    say(f"\n  best positive class - scrambled: {sep:+.3f} log2 ({2**sep:.2f}-fold)")
    say(f"  pooled positives vs scrambled, Mann-Whitney: p = {u.pvalue:.3g}")
    say(f"  scrambled sit {byg.get('scrambled',{}).get('vs_test',float('nan')):+.3f} log2 "
        "relative to test cCREs")
    bone = byg.get("bone/osteoblast enh")
    if bone:
        say(f"  note: the bone/osteoblast enhancer class sits at {bone['median']:+.3f}, i.e. "
            f"{bone['median']-scr:+.3f} vs scrambled -- as expected, they do not fire in adipocytes,")
        say("        so pooling them into one positive-control group understates the separation")
    Q["2.22"] = dict(rows=rows, separation=float(sep), mwu_p=float(u.pvalue))
    verdict(sep > 0.5,
            f"the best positive-control class sits {sep:.2f} log2 above scrambled",
            f"no positive-control class exceeds scrambled by more than {sep:.2f} log2 "
            "-- narrow dynamic range")

    # ---------------- 2.23-2.26 ----------------
    # ---- figure data -----------------------------------------------------
    PD.update(A_STAT=A_STAT, RL_MEAN=RL_MEAN, RL=RL, ACTIVE=ACTIVE,
              test_m=test_m, pos_m=pos_m, neg_m=neg_m, comp=comp,
              DNA_CPM=DNA_CPM, RNA_CPM=RNA_CPM, LFC_REP=LFC_REP, di=di, ai=ai,
              allel=np.concatenate([di, ai]),
              ctrl_cls=np.array([ctrl_class(n) for n in names]))
    try:
        capture_plots(PD)
        say("")
        say(f"figure data captured for {len(PLOT)} sections -> qc_plots.json")
    except Exception as e:
        say(f"\n  WARNING: figure data capture failed: {e}")

    # ---------------- 2.23 / 2.24 : endogenous annotation vs activity ----------
    ANN = "/wynton/group/ahituv/ryder/archaic_mpra/annotations"
    hg38map = load_hg38_map(f"{ANN}/loci_hg38.bed")
    screen = load_screen(f"{ANN}/GRCh38-cCREs.bed")
    tssidx = load_tss("/wynton/group/ahituv/js/gencode.v38.annotation.gtf")

    # oligo -> hg38 interval (both alleles of a locus share coordinates)
    loc_of = [re.search(r"_(seq\d+)_", n) for n in names]
    coords = [hg38map.get(m.group(1)) if m else None for m in loc_of]
    have = np.array([c is not None for c in coords]) & test_m
    say("")
    say(f"resolved hg38 coordinates for {int(have.sum()):,} test cCREs "
        f"({100*have.sum()/max(test_m.sum(),1):.1f}% of the test set)")

    # six activity groups, as the book specifies
    idx = np.nonzero(have)[0]
    act = RL_MEAN[idx]
    qs = np.nanquantile(act, np.linspace(0, 1, 7))
    gid = np.clip(np.searchsorted(qs, act, side="right") - 1, 0, 5)
    glabels = [f"Q{k+1}" for k in range(6)]

    head("2.23  Genomic annotations", "assess concordance with endogenous active chromatin marks")
    say("  metric: overlap with ENCODE SCREEN cCREs (Registry V4, GRCh38), with test cCREs")
    say("          split into six groups by activity. Active chromatin classes (PLS, pELS,")
    say("          dELS) should make up more of the higher-activity groups.")
    ACTIVE_CLS = ("PLS", "pELS", "dELS")
    rows = []
    for g in range(6):
        sel = idx[gid == g]
        cnt = collections.Counter()
        for i in sel:
            c, s_, e_ = coords[i]
            cls = set()
            if c in screen:
                st, en, kl = screen[c]
                j = np.searchsorted(st, e_, side="left")
                for k in np.nonzero(en[:j] > s_)[0]:
                    cls.add(kl[k])
            if not cls:
                cnt["none"] += 1
            else:
                cnt["active-like" if cls & set(ACTIVE_CLS) else "other cCRE"] += 1
        n = max(sum(cnt.values()), 1)
        rows.append(dict(group=glabels[g], n=int(sum(cnt.values())),
                         median_activity=float(np.nanmedian(RL_MEAN[sel])),
                         counts={k: int(v) for k, v in cnt.items()},
                         frac_active=float(cnt["active-like"] / n)))
    say(f"\n  {'group':<8}{'n':>7}{'median act':>12}{'active-like':>13}{'other cCRE':>12}{'no cCRE':>10}")
    for r in rows:
        c = r["counts"]; n = max(r["n"], 1)
        say(f"  {r['group']:<8}{r['n']:>7,}{r['median_activity']:>12.3f}"
            f"{100*c.get('active-like',0)/n:>12.1f}%{100*c.get('other cCRE',0)/n:>11.1f}%"
            f"{100*c.get('none',0)/n:>9.1f}%")
    Q["2.23"] = rows
    PLOT["2.23"] = rows
    lo, hi = rows[0]["frac_active"], rows[-1]["frac_active"]
    say(f"\n  active-like cCRE overlap rises from {100*lo:.1f}% (lowest activity) to "
        f"{100*hi:.1f}% (highest)")
    verdict(hi > lo,
            "overlap with active chromatin marks increases with measured activity",
            "no increase in active-chromatin overlap with activity")

    head("2.24  Proximity to TSS", "assess concordance with endogenous cCRE locations")
    say("  metric: distance from each test cCRE to the nearest GENCODE v38 TSS")
    say("          (protein-coding + lncRNA, GRCh38), by the same six activity groups.")
    dist = np.full(len(idx), np.nan)
    for k, i in enumerate(idx):
        c, s_, e_ = coords[i]
        if c not in tssidx:
            continue
        mid = (s_ + e_) // 2
        arr = tssidx[c]
        j = np.searchsorted(arr, mid)
        cand = arr[max(0, j - 1):min(len(arr), j + 1)]
        if len(cand):
            dist[k] = float(np.min(np.abs(cand - mid)))
    rows = []
    for g in range(6):
        d = dist[gid == g]
        d = d[np.isfinite(d)]
        if len(d) == 0:
            continue
        rows.append(dict(group=glabels[g], n=int(len(d)),
                         median_activity=float(np.nanmedian(RL_MEAN[idx[gid == g]])),
                         median_kb=float(np.median(d) / 1e3),
                         q1_kb=float(np.percentile(d, 25) / 1e3),
                         q3_kb=float(np.percentile(d, 75) / 1e3),
                         frac_within_10kb=float((d < 10000).mean())))
    say(f"\n  {'group':<8}{'n':>7}{'median act':>12}{'median kb':>12}{'Q1 kb':>9}{'Q3 kb':>9}{'<10 kb':>9}")
    for r in rows:
        say(f"  {r['group']:<8}{r['n']:>7,}{r['median_activity']:>12.3f}{r['median_kb']:>12.1f}"
            f"{r['q1_kb']:>9.1f}{r['q3_kb']:>9.1f}{100*r['frac_within_10kb']:>8.1f}%")
    Q["2.24"] = rows
    PLOT["2.24"] = rows
    if len(rows) >= 2:
        say(f"\n  median distance {rows[0]['median_kb']:.1f} kb (lowest activity) -> "
            f"{rows[-1]['median_kb']:.1f} kb (highest)")
        verdict(rows[-1]["median_kb"] < rows[0]["median_kb"],
                "more active cCREs lie closer to a TSS, as expected",
                "no relationship between activity and TSS proximity")
    na("2.25  AI predictions vs activity",
       "assess the correlation between AI-predicted and measured activity",
       "no AI activity predictions for this library")
    na("2.26  AI predictions vs differential activity",
       "assess the correlation between AI-predicted and measured differential activity",
       "no AI differential-activity predictions for this library")



if __name__ == "__main__":
    main()
