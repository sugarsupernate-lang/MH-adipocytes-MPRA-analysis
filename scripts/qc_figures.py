#!/usr/bin/env python3
"""Themed SVG figures for the MPRA QC report.

Reads output/qc/qc_plots.json (written by mpra_qc_suite.py) and renders one
figure per QC-book section.  Colours are CSS custom properties so a figure
reads correctly in both the light and dark artifact themes.
"""
import json
import math

W, H = 470, 262
L, R, T, B = 54, 16, 16, 44
PW, PH = W - L - R, H - T - B
FS = 'font-family="IBM Plex Mono, monospace" font-size="9"'


def esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _sx(v, lo, hi):
    return L + (v - lo) / (hi - lo) * PW if hi > lo else L


def _sy(v, lo, hi):
    return T + PH - (v - lo) / (hi - lo) * PH if hi > lo else T + PH


def _ticks(lo, hi, n=5):
    if hi <= lo:
        return [lo]
    raw = (hi - lo) / n
    mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
    step = min([1, 2, 2.5, 5, 10], key=lambda m: abs(m * mag - raw)) * mag
    t, out = math.ceil(lo / step) * step, []
    while t <= hi + step * 1e-9:
        out.append(round(t, 10))
        t += step
    return out or [lo, hi]


def _fmt(v):
    a = abs(v)
    if a >= 1e6: return f"{v/1e6:g}M"
    if a >= 1e3: return f"{v/1e3:g}k"
    if a == 0: return "0"
    if a < 0.01: return f"{v:.3f}"
    if a < 1: return f"{v:.2f}"
    return f"{v:g}"


def frame(xlo, xhi, ylo, yhi, xlab, ylab, xt=None, yt=None, xfmt=_fmt, yfmt=_fmt):
    p = [f'<rect x="{L}" y="{T}" width="{PW}" height="{PH}" fill="var(--surface-2)"/>']
    for v in (yt if yt is not None else _ticks(ylo, yhi)):
        if v < ylo - 1e-9 or v > yhi + 1e-9: continue
        y = _sy(v, ylo, yhi)
        p.append(f'<line x1="{L}" y1="{y:.1f}" x2="{L+PW}" y2="{y:.1f}" stroke="var(--line-soft)" stroke-width="1"/>')
        p.append(f'<text x="{L-6}" y="{y+3:.1f}" text-anchor="end" fill="var(--muted)" {FS}>{esc(yfmt(v))}</text>')
    for v in (xt if xt is not None else _ticks(xlo, xhi)):
        if v < xlo - 1e-9 or v > xhi + 1e-9: continue
        x = _sx(v, xlo, xhi)
        p.append(f'<text x="{x:.1f}" y="{T+PH+14}" text-anchor="middle" fill="var(--muted)" {FS}>{esc(xfmt(v))}</text>')
    p.append(f'<rect x="{L}" y="{T}" width="{PW}" height="{PH}" fill="none" stroke="var(--line)" stroke-width="1"/>')
    p.append(f'<text x="{L+PW/2:.0f}" y="{H-8}" text-anchor="middle" fill="var(--muted)" {FS}>{esc(xlab)}</text>')
    p.append(f'<text x="12" y="{T+PH/2:.0f}" text-anchor="middle" fill="var(--muted)" {FS} '
             f'transform="rotate(-90 12 {T+PH/2:.0f})">{esc(ylab)}</text>')
    return p


def wrap(parts, cls=""):
    return (f'<svg viewBox="0 0 {W} {H}" role="img" class="{cls}" '
            f'xmlns="http://www.w3.org/2000/svg">' + "".join(parts) + "</svg>")


def legend(items, y=None):
    """items: [(label, color, kind)] kind in bar|line|dot"""
    out, x = [], L + 4
    yy = T + 11 if y is None else y
    for lab, col, kind in items:
        if kind == "line":
            out.append(f'<line x1="{x}" y1="{yy-3}" x2="{x+13}" y2="{yy-3}" stroke="{col}" stroke-width="2.2"/>')
        else:
            out.append(f'<rect x="{x}" y="{yy-7}" width="13" height="8" fill="{col}" opacity="{0.85 if kind=="bar" else 1}"/>')
        out.append(f'<text x="{x+17}" y="{yy}" fill="var(--muted)" {FS}>{esc(lab)}</text>')
        x += 22 + 5.4 * len(lab)
    return out


# ---------------------------------------------------------------- primitives
def fig_hist(series, xlab, ylab="cCREs", logy=False):
    """series: [(hist1d, colour, label)]"""
    series = [(h, c, l) for h, c, l in series if h]
    if not series: return None
    xlo = min(h["edges"][0] for h, _, _ in series)
    xhi = max(h["edges"][-1] for h, _, _ in series)
    top = max(max(h["counts"]) for h, _, _ in series) or 1
    ylo, yhi = (0, top * 1.1) if not logy else (0, math.log10(top) * 1.1)
    p = frame(xlo, xhi, ylo, yhi, xlab, ylab)
    for h, col, _ in series:
        e, c = h["edges"], h["counts"]
        for i, n in enumerate(c):
            if n <= 0: continue
            x0, x1 = _sx(e[i], xlo, xhi), _sx(e[i + 1], xlo, xhi)
            v = math.log10(n) if logy else n
            y = _sy(v, ylo, yhi)
            p.append(f'<rect x="{x0:.1f}" y="{y:.1f}" width="{max(x1-x0-0.4,0.6):.1f}" '
                     f'height="{T+PH-y:.1f}" fill="{col}" opacity="0.75"/>')
    lg = [(l, c, "bar") for _, c, l in series if l]
    if lg: p += legend(lg)
    return wrap(p)


def fig_density(h2, xlab, ylab, diag=False, fit=False, extra=None):
    if not h2: return None
    xlo, xhi, ylo, yhi = h2["xlo"], h2["xhi"], h2["ylo"], h2["yhi"]
    p = frame(xlo, xhi, ylo, yhi, xlab, ylab)
    n = h2["bins"]
    cw, ch = PW / n, PH / n
    vmax = math.log10(h2["vmax"] + 1) or 1
    for i, j, c in h2["cells"]:
        o = 0.10 + 0.85 * (math.log10(c + 1) / vmax)
        x = L + i * cw
        y = T + PH - (j + 1) * ch
        p.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{cw+0.4:.2f}" height="{ch+0.4:.2f}" '
                 f'fill="var(--accent)" opacity="{o:.3f}"/>')
    if diag:
        d0, d1 = max(xlo, ylo), min(xhi, yhi)
        p.append(f'<line x1="{_sx(d0,xlo,xhi):.1f}" y1="{_sy(d0,ylo,yhi):.1f}" '
                 f'x2="{_sx(d1,xlo,xhi):.1f}" y2="{_sy(d1,ylo,yhi):.1f}" '
                 f'stroke="var(--muted)" stroke-width="1.3" stroke-dasharray="4 3"/>')
    if extra: p += extra(xlo, xhi, ylo, yhi)
    return wrap(p)


def fig_lines(series, xlab, ylab, xt=None, yt=None, pct=False, pts=True):
    """series: [(xs, ys, colour, label, dashed)]"""
    series = [s for s in series if s and s[0] and s[1]]
    if not series: return None
    xlo = min(min(s[0]) for s in series); xhi = max(max(s[0]) for s in series)
    ylo = min(min(s[1]) for s in series); yhi = max(max(s[1]) for s in series)
    pad = (yhi - ylo) * 0.12 or abs(yhi) * 0.12 or 1
    ylo, yhi = ylo - pad, yhi + pad
    f = (lambda v: f"{100*v:g}%") if pct else _fmt
    p = frame(xlo, xhi, ylo, yhi, xlab, ylab, xt=xt, yt=yt, yfmt=f)
    for xs, ys, col, lab, dash in series:
        d = " ".join(f"{_sx(x,xlo,xhi):.1f},{_sy(y,ylo,yhi):.1f}" for x, y in zip(xs, ys))
        da = ' stroke-dasharray="5 3"' if dash else ''
        p.append(f'<polyline points="{d}" fill="none" stroke="{col}" stroke-width="2.1" '
                 f'stroke-linejoin="round"{da}/>')
        if pts and len(xs) <= 20:
            for x, y in zip(xs, ys):
                p.append(f'<circle cx="{_sx(x,xlo,xhi):.1f}" cy="{_sy(y,ylo,yhi):.1f}" r="2.5" fill="{col}"/>')
    lg = [(l, c, "line") for _, _, c, l, _ in series if l]
    if lg: p += legend(lg)
    return wrap(p)


def fig_box(boxes, ylab, xlab=""):
    boxes = [b for b in boxes if b]
    if not boxes: return None
    ylo = min(b["lo"] for b in boxes); yhi = max(b["hi"] for b in boxes)
    pad = (yhi - ylo) * 0.1 or 1
    ylo, yhi = ylo - pad, yhi + pad
    p = frame(0, len(boxes), ylo, yhi, xlab, ylab, xt=[])
    bw = PW / len(boxes)
    cols = ["var(--ok)", "var(--crit)", "var(--accent)", "var(--warn)"]
    for k, b in enumerate(boxes):
        cx = L + bw * (k + 0.5)
        col = cols[k % len(cols)]
        w = min(bw * 0.42, 34)
        y1, y3 = _sy(b["q3"], ylo, yhi), _sy(b["q1"], ylo, yhi)
        p.append(f'<line x1="{cx}" y1="{_sy(b["hi"],ylo,yhi):.1f}" x2="{cx}" y2="{_sy(b["lo"],ylo,yhi):.1f}" stroke="{col}" stroke-width="1.2"/>')
        p.append(f'<rect x="{cx-w:.1f}" y="{y1:.1f}" width="{2*w:.1f}" height="{max(y3-y1,1):.1f}" fill="{col}" opacity="0.28" stroke="{col}" stroke-width="1.4"/>')
        p.append(f'<line x1="{cx-w:.1f}" y1="{_sy(b["med"],ylo,yhi):.1f}" x2="{cx+w:.1f}" y2="{_sy(b["med"],ylo,yhi):.1f}" stroke="{col}" stroke-width="2.4"/>')
        p.append(f'<text x="{cx:.0f}" y="{T+PH+15}" text-anchor="middle" fill="var(--muted)" {FS}>{esc(b["label"])}</text>')
        p.append(f'<text x="{cx:.0f}" y="{T+PH+26}" text-anchor="middle" fill="var(--muted)" font-family="IBM Plex Mono, monospace" font-size="8">n={b["n"]:,}</text>')
    return wrap(p)


def fig_binned(rows, xlab, key="per"):
    """The book's paired view: reads-per-cCRE line over a cCRE-count histogram."""
    if not rows: return None
    n = len(rows)
    top = max(r[key] for r in rows) * 1.15 or 1
    p = frame(0, n, 0, top, xlab, "DNA reads / cCRE", xt=[])
    cmax = max(r["n"] for r in rows) or 1
    bw = PW / n
    for k, r in enumerate(rows):
        h = (r["n"] / cmax) * PH * 0.42
        x = L + bw * k
        p.append(f'<rect x="{x+1:.1f}" y="{T+PH-h:.1f}" width="{bw-2:.1f}" height="{h:.1f}" '
                 f'fill="var(--muted)" opacity="0.25"/>')
    d = " ".join(f"{L+bw*(k+0.5):.1f},{_sy(r[key],0,top):.1f}" for k, r in enumerate(rows))
    p.append(f'<polyline points="{d}" fill="none" stroke="var(--accent)" stroke-width="2.2" stroke-linejoin="round"/>')
    for k, r in enumerate(rows):
        p.append(f'<circle cx="{L+bw*(k+0.5):.1f}" cy="{_sy(r[key],0,top):.1f}" r="2.4" fill="var(--accent)"/>')
    step = max(1, n // 7)
    for k, r in enumerate(rows):
        if k % step: continue
        lab = f'{r["lo"]:.2f}' if r.get("hi") is not None and r["hi"] <= 1.01 else f'{int(r["lo"])}'
        p.append(f'<text x="{L+bw*(k+0.5):.1f}" y="{T+PH+14}" text-anchor="middle" fill="var(--muted)" {FS}>{esc(lab)}</text>')
    p += legend([("reads / cCRE", "var(--accent)", "line"), ("cCREs per bin", "var(--muted)", "bar")])
    return wrap(p)


def fig_qq(d):
    if not d: return None
    p = frame(0, 1, 0, 1, "expected quantile (uniform)", "observed p-value")
    p.append(f'<line x1="{L}" y1="{T+PH}" x2="{L+PW}" y2="{T}" stroke="var(--muted)" stroke-width="1.3" stroke-dasharray="4 3"/>')
    for key, col, lab in (("test", "var(--accent)", "test cCREs"), ("neg", "var(--crit)", "negative controls")):
        s = d.get(key)
        if not s: continue
        pts = " ".join(f"{_sx(e,0,1):.1f},{_sy(o,0,1):.1f}" for e, o in zip(s["exp"], s["obs"]))
        p.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="2.1"/>')
    p += legend([("test cCREs", "var(--accent)", "line"), ("negative controls", "var(--crit)", "line"),
                 ("uniform", "var(--muted)", "line")])
    return wrap(p)


def fig_volcano(d):
    if not d: return None
    xs = [a for a, _, _ in d["pts"]]; ys = [b for _, b, _ in d["pts"]]
    if not xs: return None
    xlo, xhi = min(xs), max(xs)
    m = max(abs(xlo), abs(xhi)); xlo, xhi = -m * 1.05, m * 1.05
    yhi = max(ys) * 1.08 or 1
    p = frame(xlo, xhi, 0, yhi, "logFC   log2(Derived / Ancestral)", "-log10 FDR")
    cy = _sy(d["cut"], 0, yhi)
    p.append(f'<line x1="{L}" y1="{cy:.1f}" x2="{L+PW}" y2="{cy:.1f}" stroke="var(--muted)" stroke-width="1.2" stroke-dasharray="4 3"/>')
    p.append(f'<line x1="{_sx(0,xlo,xhi):.1f}" y1="{T}" x2="{_sx(0,xlo,xhi):.1f}" y2="{T+PH}" stroke="var(--line)" stroke-width="1"/>')
    for a, b, sig in d["pts"]:
        col = ("var(--accent)" if a > 0 else "var(--crit)") if sig else "var(--muted)"
        p.append(f'<circle cx="{_sx(a,xlo,xhi):.1f}" cy="{_sy(min(b,yhi),0,yhi):.1f}" r="1.5" '
                 f'fill="{col}" opacity="{0.75 if sig else 0.28}"/>')
    p += legend([("Derived-up", "var(--accent)", "dot"), ("Ancestral-up", "var(--crit)", "dot"),
                 ("n.s.", "var(--muted)", "dot")])
    return wrap(p)


def fig_pca(d):
    if not d: return None
    x, y = d["pc1"], d["pc2"]
    mx = max(abs(v) for v in x) * 1.25 or 1
    my = max(abs(v) for v in y) * 1.5 or 1
    p = frame(-mx, mx, -my, my, f'PC1  ({100*d["var"][0]:.1f}%)', f'PC2  ({100*d["var"][1]:.1f}%)')
    p.append(f'<line x1="{L}" y1="{_sy(0,-my,my):.1f}" x2="{L+PW}" y2="{_sy(0,-my,my):.1f}" stroke="var(--line)" stroke-width="1"/>')
    p.append(f'<line x1="{_sx(0,-mx,mx):.1f}" y1="{T}" x2="{_sx(0,-mx,mx):.1f}" y2="{T+PH}" stroke="var(--line)" stroke-width="1"/>')
    for lab, a, b in zip(d["labels"], x, y):
        col = "var(--accent)" if lab.startswith("DNA") else "var(--crit)"
        cx, cy = _sx(a, -mx, mx), _sy(b, -my, my)
        p.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="5" fill="{col}" opacity="0.85"/>')
        p.append(f'<text x="{cx:.1f}" y="{cy-9:.1f}" text-anchor="middle" fill="var(--ink)" {FS}>{esc(lab)}</text>')
    p += legend([("DNA", "var(--accent)", "dot"), ("RNA", "var(--crit)", "dot")])
    return wrap(p)


def fig_grid(rows):
    """2.19 noise grid: mean replicate r over filter x min-DNA."""
    if not rows: return None
    filts = sorted({r["filter"] for r in rows})
    mins = sorted({r["min_dna"] for r in rows})
    vals = [r["mean_r"] for r in rows]
    vlo, vhi = min(vals), max(vals)
    p = [f'<rect x="{L}" y="{T}" width="{PW}" height="{PH}" fill="var(--surface-2)"/>']
    cw, ch = PW / len(mins), PH / len(filts)
    for r in rows:
        i, j = mins.index(r["min_dna"]), filts.index(r["filter"])
        o = 0.15 + 0.8 * ((r["mean_r"] - vlo) / (vhi - vlo) if vhi > vlo else 0.5)
        x, y = L + i * cw, T + j * ch
        p.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{cw:.1f}" height="{ch:.1f}" fill="var(--accent)" opacity="{o:.3f}"/>')
        p.append(f'<text x="{x+cw/2:.1f}" y="{y+ch/2+3:.1f}" text-anchor="middle" fill="var(--ink)" '
                 f'font-family="IBM Plex Mono, monospace" font-size="9">{r["mean_r"]:.3f}</text>')
    for j, f_ in enumerate(filts):
        p.append(f'<text x="{L-6}" y="{T+ch*(j+0.5)+3:.1f}" text-anchor="end" fill="var(--muted)" {FS}>{esc(f_)}</text>')
    for i, m in enumerate(mins):
        p.append(f'<text x="{L+cw*(i+0.5):.1f}" y="{T+PH+14}" text-anchor="middle" fill="var(--muted)" {FS}>{m}</text>')
    p.append(f'<rect x="{L}" y="{T}" width="{PW}" height="{PH}" fill="none" stroke="var(--line)" stroke-width="1"/>')
    p.append(f'<text x="{L+PW/2:.0f}" y="{H-8}" text-anchor="middle" fill="var(--muted)" {FS}>minimum DNA count per cCRE</text>')
    p.append(f'<text x="12" y="{T+PH/2:.0f}" text-anchor="middle" fill="var(--muted)" {FS} '
             f'transform="rotate(-90 12 {T+PH/2:.0f})">outlier filter</text>')
    return wrap(p)


def fig_stacked(groups, classes, colors, xlab, ylab="fraction of cCREs"):
    """groups: [(label, {class: count})] rendered as stacked fractional bars."""
    if not groups: return None
    p = frame(0, len(groups), 0, 1, xlab, ylab, xt=[],
              yt=[0, .25, .5, .75, 1], yfmt=lambda v: f"{100*v:g}%")
    bw = PW / len(groups)
    for k, (lab, d) in enumerate(groups):
        tot = sum(d.get(c, 0) for c in classes) or 1
        y0 = T + PH
        x = L + bw * k + bw * 0.12
        w = bw * 0.76
        for c, col in zip(classes, colors):
            h = (d.get(c, 0) / tot) * PH
            if h <= 0: continue
            p.append(f'<rect x="{x:.1f}" y="{y0-h:.1f}" width="{w:.1f}" height="{h:.1f}" '
                     f'fill="{col}" opacity="0.85"/>')
            y0 -= h
        p.append(f'<text x="{x+w/2:.1f}" y="{T+PH+14}" text-anchor="middle" '
                 f'fill="var(--muted)" {FS}>{esc(lab)}</text>')
        p.append(f'<text x="{x+w/2:.1f}" y="{T+PH+25}" text-anchor="middle" '
                 f'fill="var(--muted)" font-family="IBM Plex Mono, monospace" font-size="8">n={tot:,}</text>')
    p += legend([(c, col, "bar") for c, col in zip(classes, colors)])
    return wrap(p)


def fig_grouped_bars(groups, series, colors, xlab, ylab, total=None):
    """groups: [(label, {series: value})]  -- grouped vertical bars."""
    if not groups: return None
    top = max(max(d.get(s, 0) for s in series) for _, d in groups) or 1
    top *= 1.18
    p = frame(0, len(groups), 0, top, xlab, ylab, xt=[])
    gw = PW / len(groups)
    bw = gw * 0.72 / len(series)
    for k, (lab, d) in enumerate(groups):
        x0 = L + gw * k + gw * 0.14
        for m, (s, col) in enumerate(zip(series, colors)):
            v = d.get(s, 0)
            h = (v / top) * PH
            x = x0 + m * bw
            p.append(f'<rect x="{x:.1f}" y="{T+PH-h:.1f}" width="{bw*0.86:.1f}" '
                     f'height="{h:.1f}" fill="{col}" opacity="0.85"/>')
            p.append(f'<text x="{x+bw*0.43:.1f}" y="{T+PH-h-3:.1f}" text-anchor="middle" '
                     f'fill="var(--ink)" font-family="IBM Plex Mono, monospace" font-size="9">{v}</text>')
        p.append(f'<text x="{L+gw*(k+0.5):.1f}" y="{T+PH+14}" text-anchor="middle" '
                 f'fill="var(--muted)" {FS}>{esc(lab)}</text>')
    if total:
        p.append(f'<text x="{L+PW-4}" y="{T+11}" text-anchor="end" fill="var(--muted)" {FS}>'
                 f'of {total} tested bases</text>')
    p += legend([(s, c, "bar") for s, c in zip(series, colors)])
    return wrap(p)
