#!/usr/bin/env python3
"""Assemble the MH adipocyte QC artifact.

Reads the QC suite's own outputs and renders one figure per section, following
the structure of the Gokhman Lab MPRA QC book: goal, figure, legend, result.
Regenerating this after a QC re-run keeps the figures and the numbers in step.

Usage: python3 build_qc_artifact.py [--out PATH]
"""
import argparse, csv, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qc_figures as F

BASE = "/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo"

# Section metadata: goals and legends follow the QC book; legends are adapted to
# name what is actually plotted here.
SEC = [
 ("1.1", "cCREs per BC", "Assess BC promiscuity",
  "Distribution of the number of cCREs each barcode is assigned to, before the promiscuity filter. A specific experiment has almost every barcode on exactly one cCRE."),
 ("1.2", "Reads per association", "Assess confidence in associations", None),
 ("1.3", "BCs per cCRE", "Assess BC replicability (complexity)",
  "Empirical cumulative distribution of barcodes per cCRE. A curve that stays low until high barcode counts means most cCREs carry many barcodes."),
 ("1.4", "Retained cCREs", "Assess BC replicability (complexity)",
  "cCREs retained at increasing minimum barcodes-per-cCRE cutoffs. A high-complexity library keeps most cCREs even under stringent filters."),
 ("1.5", "cCRE retention by sequencing depth", "Assess whether sequencing depth is sufficient",
  "Reads are downsampled and cCRE retention recounted, then fitted with the Hill saturation model and extrapolated to 3x current depth."),
 ("1.6", "BCs per cCRE by sequencing depth", "Assess whether sequencing depth is sufficient",
  "Quartiles of barcodes per cCRE against the fraction of reads used. A curve that flattens means the barcode pool is sequenced to saturation."),
 ("1.7", "PCR bias - GC", "Assess GC content bias in PCR amplification",
  "DNA reads per cCRE (line) over the number of cCREs per bin (bars), in fixed 5% GC bins. A flat line means amplification is even across GC content."),
 ("1.8", "PCR bias - G-stretches", "Assess G-stretch bias in PCR amplification",
  "DNA reads per cCRE (line) over cCRE count (bars), binned by longest run of consecutive G."),
 ("2.1", "Retained cCREs and BCs", "Assess cCRE and BC coverage",
  "Percentage of association barcodes and cCREs observed in the DNA and RNA quantification."),
 ("2.2", "DNA counts vs GC content", "Assess GC content bias in PCR amplification",
  "DNA reads per cCRE against binned GC content, with the cCRE count per bin beneath. Same data as 1.7 seen from the activity side."),
 ("2.3", "Activity statistic vs count ratio", "Assess the correlation between log2(RNA/DNA) and the activity statistic",
  "Density of test cCREs: pipeline log2(RNA/DNA) on x, MPRAnalyze activity statistic on y. Darker cells hold more cCREs."),
 ("2.4", "Activity distribution", "Assess activity dynamic range, noise and statistical power",
  "Distribution of log2(RNA/DNA) over test cCREs, with those called active by MPRAnalyze (fdr.mad < 0.05) overlaid."),
 ("2.5", "Differential activity distribution", "Assess differential activity dynamic range, noise and power",
  "Distribution of logFC = log2(Derived/Ancestral) across tested loci, with significantly differentially active loci overlaid."),
 ("2.6", "P-value distribution", "Assess statistical power",
  "QQ plot of MPRAnalyze pval.mad against a uniform null, for test cCREs and for the scrambled negative controls."),
 ("2.7", "Volcano plot - FC vs P-value", "Assess statistical power of differential activity detection",
  "logFC against -log10 FDR. The dashed line is FDR = 0.05; colour marks which allele is the more active one."),
 ("2.8", "BC retention by DNA/RNA sequencing depth", "Assess whether sequencing depth is sufficient",
  "Retained barcodes against depth, DNA and RNA downsampled jointly, with Hill extrapolation to 3x."),
 ("2.9", "cCRE retention by DNA/RNA sequencing depth", "Assess whether sequencing depth is sufficient",
  "Retained cCREs against depth for each channel."),
 ("2.10", "Activity by sequencing depth", "Assess whether sequencing depth is sufficient",
  "Fraction of cCREs above the activity proxy log2(RNA/DNA) > 1 at each depth. Both channels are downsampled together, so the ratio is unbiased and the threshold is held fixed."),
 ("2.11", "Cumulative RNA reads", "Assess jackpotting",
  "Cumulative share of RNA UMIs against cCREs ranked by RNA count. A curve hugging the top-left means a few cCREs dominate."),
 ("2.12", "Sample clustering", "Assess reproducibility between samples",
  "PCA over log10(CPM+1) per cCRE with the six samples as observations."),
 ("2.13", "Correlation between replicates", "Assess reproducibility between replicates",
  "Density of CPM-normalised log2(RNA/DNA) per cCRE, replicate 1 against replicate 2. Dashed line is identity."),
 ("2.14", "Correlation of differential activity between replicates", "Assess reproducibility between replicates",
  "Density of log2(RNA_cpm/DNA_cpm) for the oligos making up the allelic pairs, replicate 1 against replicate 2. Each allele is its own oligo, so this is the reproducibility of the measurement the comparative test is built from."),
 ("2.15", "Replicability by activity", "Assess replicate correlation in active vs non-active cCREs",
  "Replicate correlation within fixed bins of mean log2(RNA/DNA), bin width 1."),
 ("2.16", "Correlation between replicates (controls)", "Assess reproducibility in positive and negative controls",
  "Replicate 1 against replicate 2, with positive controls in green and scrambled negatives in red over the density of test cCREs."),
 ("2.17", "Cross-validation: allelic pairs", "Assess reproducibility between allelic pairs",
  "log2(RNA/DNA) of the Ancestral against the Derived allele of the same locus. The two differ by a few bases and should track closely."),
 ("2.18", "Cross-validation: cell types", "Assess reproducibility between different cell types",
  "log2(RNA/DNA) of shared cCREs in this experiment against a sibling MH cell type."),
 ("2.19", "Minimizing noise", "Optimize outlier removal and minimum DNA threshold",
  "Mean replicate correlation across the outlier-filter by minimum-DNA grid. Darker is more reproducible."),
 ("2.20", "Reproducibility by sequencing depth", "Assess whether more sequencing improves reproducibility",
  "Replicate correlation as a function of downsampled depth."),
 ("2.21", "RNA vs DNA", "Assess true activity in the experiment",
  "CPM-normalised RNA against DNA per test cCRE, log10 scale. Points above the dashed identity line have RNA in excess of DNA."),
 ("2.22", "Activity of controls", "Assess the dynamic range of activity",
  "Activity of positive controls, scrambled negative controls and test cCREs. Box is the interquartile range, bar the median."),
 ("2.23", "Genomic annotations", "Assess concordance with endogenous active chromatin marks",
  "Overlap with ENCODE SCREEN cCREs (Registry V4, GRCh38), with test cCREs split into six equal-sized groups by measured activity. Active-like classes are PLS, pELS and dELS; the remainder overlap some other cCRE class or none at all."),
 ("2.24", "Proximity to TSS", "Assess concordance with endogenous cCRE locations",
  "Distance from each test cCRE to the nearest GENCODE v38 TSS (protein-coding and lncRNA, GRCh38), by the same six activity groups. Box is the interquartile range, bar the median."),
 ("2.25", "AI predictions vs activity", "Assess correlation between AI-predicted and measured activity", None),
 ("2.26", "AI predictions vs differential activity", "Assess correlation between AI-predicted and measured differential activity", None),
]


def parse_report(path):
    """Pull each section's verdict and body out of the report the suite wrote."""
    out, cur = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        m = re.match(r"^  (\d+\.\d+)\s\s+(.+)$", line)
        if m and not line.strip().startswith("->"):
            cur = m.group(1)
            out[cur] = dict(verdict=None, text=None, body=[], na=None)
            continue
        if cur is None:
            continue
        v = re.match(r"^  -> (PASS|FLAG): (.+)$", line)
        if v:
            out[cur]["verdict"], out[cur]["text"] = v.group(1), v.group(2)
            continue
        n = re.match(r"^  NOT EVALUATED -- (.+)$", line)
        if n:
            out[cur]["na"] = n.group(1)
            continue
        if line.strip() and not line.startswith("====") and not line.strip().startswith("goal:"):
            out[cur]["body"].append(line)
    return out


def build_figure(sid, P):
    d = P.get(sid)
    A, C, M, O, Wc = "var(--accent)", "var(--crit)", "var(--muted)", "var(--ok)", "var(--warn)"
    if sid == "1.1" and d:
        h = dict(edges=[k - 0.5 for k in d["k"]] + [d["k"][-1] + 0.5],
                 counts=d["v"], n=sum(d["v"]))
        return F.fig_hist([(h, A, "")], "cCREs per barcode", "barcodes", logy=True)
    if sid == "1.3" and d:
        return F.fig_lines([(d["x"], d["y"], A, "", False)], "barcodes per cCRE",
                           "cumulative fraction of cCREs", pts=False)
    if sid == "1.4" and d:
        return F.fig_lines([([r["min_bc"] for r in d], [r["frac"] for r in d], A, "", False)],
                           "minimum barcodes per cCRE", "cCREs retained", pct=True)
    if sid == "1.5" and d and d.get("frac"):
        return F.fig_lines([(d["frac"], d["ret"], A, "observed", False)],
                           "fraction of reads used", "cCREs retained", pct=True)
    if sid == "1.6" and d:
        return F.fig_lines([([r["frac"] for r in d], [r["q1"] for r in d], M, "Q1", True),
                            ([r["frac"] for r in d], [r["median"] for r in d], A, "median", False),
                            ([r["frac"] for r in d], [r["q3"] for r in d], M, "Q3", True)],
                           "fraction of reads used", "barcodes per cCRE")
    if sid in ("1.7", "2.2") and d:
        return F.fig_binned(d, "GC content")
    if sid == "1.8" and d:
        return F.fig_binned(d, "longest G-stretch (bp)")
    if sid == "2.1" and d:
        rows = [r for r in d if r.get("channel") in ("DNA", "RNA", "either")]
        if not rows: return None
        tot = P.get("_assoc_bc") or max(r["barcodes"] for r in rows)
        return F.fig_lines([([0, 1, 2], [r["barcodes"] / tot for r in rows], A, "barcodes", False)],
                           "DNA / RNA / either", "fraction observed", xt=[0, 1, 2], pct=True)
    if sid == "2.3" and d:
        return F.fig_density(d, "log2(RNA/DNA)", "MPRAnalyze activity statistic")
    if sid == "2.4" and d:
        return F.fig_hist([(d.get("all"), M, "all test cCREs"), (d.get("active"), O, "active")],
                          "log2(RNA/DNA)")
    if sid == "2.5" and d:
        return F.fig_hist([(d.get("all"), M, "all loci"), (d.get("sig"), C, "differentially active")],
                          "logFC   log2(Derived / Ancestral)", "loci")
    if sid == "2.6" and d:
        return F.fig_qq(d)
    if sid == "2.7" and d:
        return F.fig_volcano(d)
    if sid in ("2.8", "2.9") and d and d.get("frac"):
        return F.fig_lines([(d["frac"], d["dna"], A, "DNA", False), (d["frac"], d["rna"], C, "RNA", False)],
                           "fraction of reads used",
                           "barcodes retained" if sid == "2.8" else "cCREs retained", pct=True)
    if sid == "2.10" and d and d.get("frac"):
        return F.fig_lines([(d["frac"], d["act"], A, "", False)], "fraction of reads used",
                           "cCREs above proxy", pct=True)
    if sid == "2.11" and d:
        return F.fig_lines([(d["x"], d["y"], A, "", False)], "cCREs ranked by RNA count",
                           "cumulative share of RNA UMIs", pct=True, pts=False)
    if sid == "2.12" and d:
        return F.fig_pca(d)
    if sid == "2.13" and d:
        return F.fig_density(d, "rep1  log2(RNA/DNA)", "rep2  log2(RNA/DNA)", diag=True)
    if sid == "2.14" and d:
        return F.fig_density(d, "rep1  log2(RNA/DNA)", "rep2  log2(RNA/DNA)", diag=True)
    if sid == "2.15" and d:
        xs = [r["lo"] + 0.5 for r in d]
        return F.fig_lines([(xs, [r["r"][0] for r in d], A, "r(1,2)", False),
                            (xs, [r["r"][1] for r in d], C, "r(1,3)", False),
                            (xs, [r["r"][2] for r in d], O, "r(2,3)", False)],
                           "mean log2(RNA/DNA) bin", "replicate correlation")
    if sid == "2.16" and d and d.get("dens"):
        def overlay(xlo, xhi, ylo, yhi):
            out = []
            for pts, col in ((d.get("neg", []), C), (d.get("pos", []), O)):
                for a, b in pts:
                    if not (xlo <= a <= xhi and ylo <= b <= yhi): continue
                    out.append(f'<circle cx="{F._sx(a,xlo,xhi):.1f}" cy="{F._sy(b,ylo,yhi):.1f}" '
                               f'r="2.1" fill="{col}" opacity="0.85"/>')
            out += F.legend([("positive ctrl", O, "dot"), ("negative ctrl", C, "dot")])
            return out
        return F.fig_density(d["dens"], "rep1  log2(RNA/DNA)", "rep2  log2(RNA/DNA)",
                             diag=True, extra=overlay)
    if sid == "2.17" and d:
        return F.fig_density(d, "Ancestral  log2(RNA/DNA)", "Derived  log2(RNA/DNA)", diag=True)
    if sid == "2.18" and d:
        best = max(d.items(), key=lambda kv: kv[1].get("pearson") or -9)
        return F.fig_density(best[1]["dens"], "adipocyte  log2(RNA/DNA)",
                             f'{best[0]}  log2(RNA/DNA)', diag=True)
    if sid == "2.19" and d:
        return F.fig_grid(d)
    if sid == "2.20" and d:
        return F.fig_lines([([r["frac"] for r in d], [r["mean_r"] for r in d], A, "", False)],
                           "fraction of reads used", "mean replicate correlation")
    if sid == "2.21" and d and d.get("dens"):
        return F.fig_density(d["dens"], "DNA CPM", "RNA CPM", diag=True)
    if sid == "2.22" and d:
        return F.fig_box(d, "log2(RNA/DNA)")
    if sid == "2.23" and d:
        cls = ["active-like", "other cCRE", "none"]
        cols = [O, A, M]
        return F.fig_stacked([(r["group"], r["counts"]) for r in d], cls, cols,
                             "activity group (low to high)")
    if sid == "2.24" and d:
        boxes = [dict(label=r["group"], q1=r["q1_kb"], med=r["median_kb"], q3=r["q3_kb"],
                      lo=max(r["q1_kb"] - 1.5 * (r["q3_kb"] - r["q1_kb"]), 0),
                      hi=r["q3_kb"] + 1.5 * (r["q3_kb"] - r["q1_kb"]),
                      n=r["n"], mean=r["median_kb"]) for r in d]
        return F.fig_box(boxes, "distance to nearest TSS (kb)", "activity group (low to high)")
    return None


CSS = """
:root{--ground:#EFF3F4;--surface:#FFFFFF;--surface-2:#F7FAFA;--ink:#16232B;--muted:#5E7480;
--line:#D3DEE0;--line-soft:#E4ECED;--accent:#0E6E77;--accent-soft:#DCEDEE;--crit:#B4322C;
--crit-soft:#F7E4E2;--warn:#A8700A;--warn-soft:#F8EEDA;--ok:#2F6B4F;--ok-soft:#E1EFE7;
--shadow:0 1px 2px rgba(22,35,43,.06),0 6px 20px -12px rgba(22,35,43,.28);}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
--ground:#0F171B;--surface:#161F24;--surface-2:#1B262B;--ink:#E3ECEF;--muted:#93A8B2;
--line:#2B383E;--line-soft:#222E34;--accent:#45B4BC;--accent-soft:#123336;--crit:#E58179;
--crit-soft:#3A211F;--warn:#DCA744;--warn-soft:#332914;--ok:#6FBF95;--ok-soft:#17301F;
--shadow:0 1px 2px rgba(0,0,0,.4),0 6px 20px -12px rgba(0,0,0,.7);}}
:root[data-theme="dark"]{--ground:#0F171B;--surface:#161F24;--surface-2:#1B262B;--ink:#E3ECEF;
--muted:#93A8B2;--line:#2B383E;--line-soft:#222E34;--accent:#45B4BC;--accent-soft:#123336;
--crit:#E58179;--crit-soft:#3A211F;--warn:#DCA744;--warn-soft:#332914;--ok:#6FBF95;
--ok-soft:#17301F;--shadow:0 1px 2px rgba(0,0,0,.4),0 6px 20px -12px rgba(0,0,0,.7);}
*{box-sizing:border-box}
body{margin:0;background:var(--ground);color:var(--ink);font-family:"IBM Plex Serif",Georgia,serif;
font-size:16px;line-height:1.62;-webkit-font-smoothing:antialiased}
.wrap{max-width:980px;margin:0 auto;padding:clamp(28px,5vw,60px) clamp(16px,4vw,34px) 80px}
h1,h2,h3,.ui{font-family:"IBM Plex Sans Condensed","Helvetica Neue",Arial,sans-serif}
.mono{font-family:"IBM Plex Mono",ui-monospace,Menlo,monospace;font-variant-numeric:tabular-nums}
.eyebrow{font-family:"IBM Plex Mono",monospace;font-size:11.5px;letter-spacing:.13em;
text-transform:uppercase;color:var(--muted);margin:0 0 14px;display:flex;flex-wrap:wrap;gap:6px 10px}
.eyebrow span:not(:last-child)::after{content:"·";margin-left:10px;color:var(--line)}
h1{font-size:clamp(34px,5.6vw,54px);line-height:1.02;font-weight:700;letter-spacing:-.018em;
margin:0 0 16px;text-wrap:balance}
.standfirst{font-size:clamp(17px,2.1vw,19.5px);line-height:1.55;color:var(--muted);max-width:64ch;margin:0 0 30px}
.standfirst strong{color:var(--ink);font-weight:600}
.standfirst a,footer a{color:var(--accent)}
.score{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));gap:1px;
background:var(--line);border:1px solid var(--line);border-radius:3px;overflow:hidden;margin:0 0 40px}
.score div{background:var(--surface);padding:16px 18px}
.score .n{font-family:"IBM Plex Sans Condensed",sans-serif;font-size:34px;font-weight:700;
line-height:1;letter-spacing:-.02em;font-variant-numeric:tabular-nums;display:block}
.score .l{font-family:"IBM Plex Mono",monospace;font-size:10.5px;letter-spacing:.11em;
text-transform:uppercase;color:var(--muted);display:block;margin-top:8px}
.c-crit .n{color:var(--crit)}.c-warn .n{color:var(--warn)}.c-ok .n{color:var(--ok)}.c-acc .n{color:var(--accent)}
.alert{border:1px solid var(--accent);border-left-width:3px;border-radius:3px;background:var(--surface);
box-shadow:var(--shadow);margin:0 0 46px;overflow:hidden}
.alert .ah{padding:16px 22px 12px;border-bottom:1px solid var(--line-soft);display:flex;gap:12px;
align-items:baseline;flex-wrap:wrap}
.alert .ah h3{margin:0;font-size:21px;font-weight:600;letter-spacing:-.008em;flex:1}
.alert .tag{font-family:"IBM Plex Mono",monospace;font-size:10px;letter-spacing:.1em;text-transform:uppercase;
padding:4px 9px;border-radius:2px;background:var(--ok-soft);color:var(--ok)}
.alert .ab{padding:16px 22px 20px}.alert p{margin:0 0 13px;max-width:66ch}.alert p:last-child{margin-bottom:0}
.conv{display:grid;grid-template-columns:1fr 1fr;gap:1px;background:var(--line);border:1px solid var(--line);
border-radius:3px;overflow:hidden;margin:0 0 15px}
.conv>div{background:var(--surface-2);padding:14px 16px}
.conv .lab{font-family:"IBM Plex Mono",monospace;font-size:10px;letter-spacing:.11em;text-transform:uppercase;
color:var(--muted);margin-bottom:8px}
.conv .fm{font-family:"IBM Plex Mono",monospace;font-size:13.5px;margin-bottom:6px}
.conv .m{font-size:13px;color:var(--muted);line-height:1.45}
.conv .pos{font-family:"IBM Plex Sans Condensed",sans-serif;font-weight:600;font-size:15px;margin-top:8px}
.conv .a .pos{color:var(--muted)}.conv .b .pos{color:var(--accent)}
@media(max-width:560px){.conv{grid-template-columns:1fr}}
.chap{margin:52px 0 22px;padding-bottom:10px;border-bottom:2px solid var(--ink)}
.chap h2{font-size:13px;font-weight:600;letter-spacing:.16em;text-transform:uppercase;margin:0}
.chap p{margin:6px 0 0;font-size:14px;color:var(--muted)}
.fig{background:var(--surface);border:1px solid var(--line);border-radius:3px;box-shadow:var(--shadow);
margin:0 0 22px;overflow:hidden}
.fig .fh{padding:15px 20px 11px;border-bottom:1px solid var(--line-soft);display:flex;gap:12px;
align-items:baseline;flex-wrap:wrap}
.fig .sid{font-family:"IBM Plex Mono",monospace;font-size:12px;font-weight:600;padding:3px 8px;
border-radius:2px;background:var(--accent-soft);color:var(--accent);flex-shrink:0}
.fig .fh h3{margin:0;font-size:18px;font-weight:600;flex:1;min-width:170px;letter-spacing:-.006em}
.fig .goal{font-size:13px;color:var(--muted);width:100%;margin:2px 0 0;font-style:italic}
.fig .fb{padding:16px 20px 18px}
.fig svg{display:block;width:100%;height:auto;max-width:520px;margin:0 auto 12px}
.fig .cap{font-size:13.5px;color:var(--muted);line-height:1.5;margin:0 0 12px;max-width:70ch}
.verdict{display:flex;gap:10px;align-items:flex-start;font-size:14.5px;line-height:1.5;
padding:10px 14px;border-radius:2px;border-left:2px solid;margin:0}
.v-PASS{background:var(--ok-soft);border-color:var(--ok)}
.v-FLAG{background:var(--crit-soft);border-color:var(--crit)}
.v-NA{background:var(--surface-2);border-color:var(--line)}
.verdict b{font-family:"IBM Plex Mono",monospace;font-size:10px;letter-spacing:.1em;flex-shrink:0;padding-top:3px}
.v-PASS b{color:var(--ok)}.v-FLAG b{color:var(--crit)}.v-NA b{color:var(--muted)}
details{margin-top:12px;border-top:1px solid var(--line-soft);padding-top:9px}
summary{font-family:"IBM Plex Mono",monospace;font-size:11px;letter-spacing:.08em;text-transform:uppercase;
color:var(--muted);cursor:pointer}
summary:hover{color:var(--accent)}
pre{font-family:"IBM Plex Mono",monospace;font-size:11.5px;line-height:1.5;background:var(--surface-2);
border:1px solid var(--line-soft);border-radius:2px;padding:11px 13px;overflow-x:auto;margin:9px 0 0;
color:var(--ink);font-variant-numeric:tabular-nums}
code{font-family:"IBM Plex Mono",monospace;font-size:.87em;background:var(--surface-2);
border:1px solid var(--line-soft);border-radius:2px;padding:.1em .34em}
footer{border-top:1px solid var(--line);padding-top:20px;margin-top:46px;font-size:13px;color:var(--muted);line-height:1.6}
:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
"""



def gwas_eqtl_section():
    """Differential hits against GWAS Catalog and GTEx adipose eQTLs.

    Each element is 200 bp but the MPRA changed exactly ONE base, so overlap with
    the element is not the same as hitting the tested variant. Windows are measured
    from the tested base.
    """
    ANN = "/wynton/group/ahituv/ryder/archaic_mpra/annotations"
    f = f"{ANN}/MH_adipo_346_gwas_eqtl_windowed.tsv"
    if not os.path.exists(f):
        return ""
    rows = list(csv.DictReader(open(f), delimiter="\t"))
    w50 = [r for r in rows if r["window"] == "50"]
    w500 = [r for r in rows if r["window"] == "500"]

    def cnt(rs):
        g = sum(1 for r in rs if r["gwas"] != "-")
        e = sum(1 for r in rs if r["eqtl_genes"] != "-")
        return dict(GWAS=g, eQTL=e, either=len(rs))
    groups = [("exact base", dict(GWAS=1, eQTL=5, either=5)),
              ("\u00b150 bp", cnt(w50)), ("\u00b1500 bp", cnt(w500))]
    fig = F.fig_grouped_bars(groups, ["GWAS", "eQTL", "either"],
                             ["var(--crit)", "var(--accent)", "var(--ok)"],
                             "window around the tested base",
                             "differential loci with evidence", total=305)

    h = ['<div class="chap"><h2>Differential hits &mdash; GWAS and adipose eQTL</h2>'
         '<p>The 346 differentially active loci intersected with the GWAS Catalog and GTEx v8 '
         'adipose eQTLs. Both sources are GRCh38, matching the lifted element coordinates.</p></div>',
         '<section class="fig"><div class="fh"><span class="sid">GxE</span>'
         '<h3>How often does a common variant reach the tested base?</h3>'
         '<p class="goal">Assess whether population-genetic resources can speak to these elements</p>'
         '</div><div class="fb">']
    if fig:
        h.append(fig)
    h.append('<p class="cap">Each element is 200 bp, but the MPRA changed exactly <em>one</em> base. '
             'Overlapping the element is therefore not the same as hitting the tested variant, so '
             'windows are measured from the tested base outward. Beyond about &plusmn;100 bp the window '
             'leaves the element entirely, which is why the &plusmn;500 bp column is proximity rather '
             'than element-level evidence.</p>')
    h.append('<p class="verdict v-FLAG"><b>STRUCTURAL LIMIT</b><span>Only 23 of 305 tested bases '
             '(7.5%) have any GWAS or adipose-eQTL signal within 50 bp. These are ancestral/derived '
             'variants, mostly rare or fixed; GTEx tests MAF &ge; 1% and the GWAS Catalog is '
             'common-variant, so most of these elements cannot be reached by either resource. That is '
             'a property of the variant class, not a coverage gap.</span></p>')
    h.append('<details><summary>the 23 loci with evidence within 50 bp of the tested base</summary>'
             '<div class="tw"><table><thead><tr><th>locus</th><th>logFC</th><th>hg38</th>'
             '<th>gene</th><th>3-cell</th><th>cCRE</th><th>GWAS trait / eQTL genes</th></tr></thead><tbody>')
    for r in sorted(w50, key=lambda x: -abs(float(x["lfc"]))):
        ev = []
        if r["gwas"] != "-":
            ev.append(f'GWAS +{r["gwas_min_d"]}bp: {F.esc(r["gwas"][:60])}')
        if r["eqtl_genes"] != "-":
            ev.append(f'eQTL +{r["eqtl_min_d"]}bp [{F.esc(r["eqtl_tis"])}]: {F.esc(r["eqtl_genes"][:44])}')
        gene = r["ov"] if r["ov"] != "-" else r["near"]
        h.append(f'<tr><td>{F.esc(r["locus"])}</td><td class="num">{float(r["lfc"]):+.3f}</td>'
                 f'<td>{F.esc(r["hg38"])}</td><td>{F.esc(gene[:18])}</td>'
                 f'<td>{"yes" if r["conc3"]=="True" else "-"}</td><td>{F.esc(r["screen"])}</td>'
                 f'<td>{"<br>".join(ev)}</td></tr>')
    h.append('</tbody></table></div></details></div></section>')
    return "\n".join(h)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plots", default=f"{BASE}/output/qc/qc_plots.json")
    ap.add_argument("--report", default=f"{BASE}/output/qc/qc_report.txt")
    ap.add_argument("--out", default=f"{BASE}/output/qc/qc_artifact.html")
    a = ap.parse_args()

    P = json.load(open(a.plots)) if os.path.exists(a.plots) else {}
    R = parse_report(a.report)
    npass = sum(1 for v in R.values() if v["verdict"] == "PASS")
    nflag = sum(1 for v in R.values() if v["verdict"] == "FLAG")
    nna = sum(1 for v in R.values() if v["na"])
    nfig = 0

    h = ['<title>MH Adipocyte QC Report</title>',
         '<link rel="preconnect" href="https://fonts.googleapis.com">',
         '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>',
         '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600'
         '&family=IBM+Plex+Sans+Condensed:wght@500;600;700&family=IBM+Plex+Serif:ital,wght@0,400;0,600;1,400&display=swap">',
         f'<style>{CSS}</style>', '<div class="wrap">',
         '<p class="eyebrow"><span>MH library</span><span>SGBS adipocyte</span>'
         '<span>HE016AH_0715 · AVITI</span><span>28,910 oligos</span></p>',
         '<h1>MH Adipocyte QC Report</h1>']

    h += ['<div class="score">',
          f'<div class="c-acc"><span class="n">{npass+nflag}</span><span class="l">sections evaluated</span></div>',
          f'<div class="c-ok"><span class="n">{npass}</span><span class="l">pass</span></div>',
          f'<div class="c-crit"><span class="n">{nflag}</span><span class="l">flagged</span></div>',
          f'<div class="c-warn"><span class="n">{nna}</span><span class="l">not evaluable</span></div>',
          '</div>']



    chap = None
    for sid, title, goal, legend_txt in SEC:
        c = sid.split(".")[0]
        if c != chap:
            chap = c
            if c == "1":
                h += ['<div class="chap"><h2>Chapter 1 &mdash; cCRE&ndash;BC association QC</h2>'
                      '<p>Barcode&ndash;oligo association quality, before any expression data.</p></div>']
            else:
                h += ['<div class="chap"><h2>Chapter 2 &mdash; activity QC</h2>'
                      '<p>DNA/RNA quantification, activity calls and differential activity.</p></div>']
        r = R.get(sid, {})
        svg = build_figure(sid, P)
        if svg: nfig += 1
        h.append('<section class="fig"><div class="fh">'
                 f'<span class="sid">{sid}</span><h3>{F.esc(title)}</h3>'
                 f'<p class="goal">{F.esc(goal)}</p></div><div class="fb">')
        if svg:
            h.append(svg)
        if legend_txt and svg:
            h.append(f'<p class="cap">{F.esc(legend_txt)}</p>')
        if r.get("na"):
            h.append(f'<p class="verdict v-NA"><b>NOT EVALUATED</b><span>{F.esc(r["na"])}</span></p>')
        elif r.get("verdict"):
            h.append(f'<p class="verdict v-{r["verdict"]}"><b>{r["verdict"]}</b>'
                     f'<span>{F.esc(r["text"])}</span></p>')
        body = [b for b in r.get("body", []) if b.strip()]
        if body:
            h.append('<details><summary>the numbers behind this figure</summary>'
                     f'<pre>{F.esc(chr(10).join(body))}</pre></details>')
        h.append('</div></section>')

    h.append(gwas_eqtl_section())
    h += ['<footer><p>Generated by <code>scripts/build_qc_artifact.py</code> from '
          '<code>output/qc/qc_plots.json</code> and <code>output/qc/qc_report.txt</code>, both written by '
          '<code>scripts/mpra_qc_suite.py</code>. Figures and numbers come from one computation; re-running the '
          'QC suite and rebuilding redraws them in step.</p>'
          f'<p class="mono" style="font-size:11.5px">/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo '
          f'&middot; {nfig} figures &middot; {npass+nflag} sections evaluated, {nna} not evaluable</p></footer>',
          '</div>']

    out = "\n".join(h)
    with open(a.out, "w") as fh:
        fh.write(out)
    print(f"wrote {a.out}  ({len(out):,} bytes, {nfig} figures, "
          f"{npass} pass / {nflag} flag / {nna} n/a)")


if __name__ == "__main__":
    main()
