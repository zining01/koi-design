#!/usr/bin/env python3
"""
koi_grid.py - compute a folding grid for a base with a grafted scale pattern.

The paper along each axis is a sequence of segments. Two modes:

PROPORTIONAL mode (every scale segment has a "length"):
  {"name": "head", "length": 0.30}                  a scale-free chunk
  {"name": "body", "length": 0.50, "scales": 22}    a chunk holding 22 scales
  {"name": "tail", "length": 0.20}
  Lengths are relative proportions (they need not sum to 1).
  The scale width follows from body length / scale count. Set either
  "paper_size_mm" (scale width is reported) or "scale_width_mm" (the paper
  size needed is reported).

  Use "folded_length" instead of "length" on a scale segment to fix the
  body's length *after* the scale pleats are collapsed (i.e. its length in
  the base without scales) instead of its length on the paper:
  {"name": "body", "folded_length": 0.50, "scales": 22}
  This needs "folded_units_per_scale": how many grid units of length one
  scale adds to the folded body (of the sum(scale_pattern) units it uses on
  the paper; the rest is hidden in the pleats). Ratios and errors are then
  measured on the folded lengths.

GRAFT mode (no scale segment has a "length"):
  {"base": L}     a piece of the original base, L in *base units*
  {"scales": n}   n scales inserted at that point
  The body length then depends on "base_units_per_scale" (how long one scale
  should be, in base units).

In both modes each scale occupies `scale_pattern` grid units (default [1, 1]:
e.g. one pleat unit + one visible unit), possibly multiplied by a whole number
so the proportions can be matched on a finer grid.

The script finds uniform N x N grids on which every chunk boundary and every
scale line falls on a grid line, then writes an SVG of the grid, a CSV of every
line position and a text summary (including how to divide the paper into N).

Usage:
  python3 koi_grid.py config.json                 # list candidate grids
  python3 koi_grid.py config.json --pick 1        # write files for candidate #1
  python3 koi_grid.py config.json --n 64          # candidate with N = 64
  python3 koi_grid.py config.json --exact         # no snapping; exact positions
"""
import argparse
import csv
import itertools
import json
import math
import os
import sys


# ---------------------------------------------------------------- config ---

def seg_length(seg):
    """Target proportional length of a segment (proportional mode)."""
    for key in ("folded_length", "length", "base"):
        if key in seg:
            return seg[key]


def folded_counts(seq, counts, cfg):
    """Segment lengths after the scale pleats of every 'folded_length'
    segment are collapsed. Other segments are measured on the paper."""
    P = sum(cfg["scale_pattern"])
    F = cfg.get("folded_units_per_scale")
    return [c * F / P if "folded_length" in s else c
            for s, c in zip(seq, counts)]


def any_folded(cfg):
    return any("folded_length" in s for a in ("x", "y") for s in cfg[a])


def seg_name(seg, i):
    if "name" in seg:
        return seg["name"]
    return ("scales%d" % i) if "scales" in seg else ("base%d" % i)


def load_config(path):
    with open(path) as f:
        cfg = json.load(f)
    cfg.setdefault("scale_pattern", [1, 1])
    cfg.setdefault("max_n", 160)
    cfg.setdefault("min_n", 8)
    if "y" not in cfg:
        cfg["y"] = cfg["x"]
    if "paper_size_mm" in cfg and "scale_width_mm" in cfg:
        sys.exit("set either 'paper_size_mm' or 'scale_width_mm', not both")
    if "paper_size_mm" not in cfg and "scale_width_mm" not in cfg:
        cfg["paper_size_mm"] = 500

    scale_segs = []
    for axis in ("x", "y"):
        if not any("scales" in s for s in cfg[axis]):
            sys.exit("axis '%s' has no scale segment" % axis)
        for seg in cfg[axis]:
            if "scales" in seg:
                if "base" in seg:
                    sys.exit("segment has both 'base' and 'scales': %r" % seg)
                if seg["scales"] < 1:
                    sys.exit("scale count must be >= 1: %r" % seg)
                scale_segs.append(seg)
            elif "folded_length" in seg:
                sys.exit("'folded_length' only applies to scale segments; "
                         "use 'length' here: %r" % seg)
            elif "base" not in seg and "length" not in seg:
                sys.exit("segment needs 'base', 'length' or 'scales': %r"
                         % seg)
    for s in scale_segs:
        if "length" in s and "folded_length" in s:
            sys.exit("give a scale segment 'length' or 'folded_length', "
                     "not both: %r" % s)
    with_len = [s for s in scale_segs
                if "length" in s or "folded_length" in s]
    if with_len and len(with_len) != len(scale_segs):
        sys.exit("either every scale segment has a 'length'/'folded_length' "
                 "(proportional mode) or none does (graft mode)")
    cfg["mode"] = "proportional" if with_len else "graft"
    if any_folded(cfg):
        F = cfg.get("folded_units_per_scale")
        P = sum(cfg["scale_pattern"])
        if F is None:
            sys.exit("'folded_length' needs 'folded_units_per_scale': how "
                     "many grid units of length one scale adds to the folded "
                     "body (out of the %g units it uses on the paper)" % P)
        if not 0 < F <= P:
            sys.exit("folded_units_per_scale must be > 0 and <= %g "
                     "(sum of scale_pattern)" % P)
    if cfg["mode"] == "graft":
        if "base_units_per_scale" not in cfg:
            sys.exit("graft mode needs 'base_units_per_scale': how long one "
                     "scale should be, in base units")
        for axis in ("x", "y"):
            for seg in cfg[axis]:
                if "length" in seg and "scales" not in seg:
                    sys.exit("graft mode: use 'base' (not 'length') for "
                             "base segments: %r" % seg)
    return cfg


# ---------------------------------------------------------------- layout ---

def layout(seq, counts, pattern):
    """
    Place every line on one axis. counts[i] = length of segment i in grid
    units (ints when snapped, floats for exact layouts). Scale segments are
    divided evenly into scales, each scale split according to `pattern`.
    Returns (lines, spans): lines = [(pos, kind)] with kind in
    edge/boundary/scale, spans = [(start, end, seg)].
    """
    P = sum(pattern)
    pos = 0.0
    lines = [(0.0, "edge")]
    spans = []
    for seg, c in zip(seq, counts):
        start = pos
        if "scales" in seg:
            unit = c / (seg["scales"] * P)
            p = start
            for _ in range(seg["scales"]):
                for u in pattern:
                    p += u * unit
                    lines.append((p, "scale"))
            lines.pop()             # the last one is the segment boundary
        pos = start + c
        lines.append((pos, "boundary"))
        spans.append((start, pos, seg))
    lines[-1] = (pos, "edge")
    return lines, spans


def scale_pitches(seq, counts, pattern):
    """Grid units per scale for every scale segment on an axis."""
    P = sum(pattern)
    return [c / s["scales"] for s, c in zip(seq, counts) if "scales" in s]


# ---------------------------------------------------- graft-mode search ---

def graft_counts(seq, k, P):
    return [s["scales"] * P if "scales" in s else int(round(k * s["base"]))
            for s in seq]


def graft_bases(seq):
    return [s["base"] for s in seq if "base" in s]


def shape_error(targets, actual):
    """Largest difference between target and actual fractions."""
    tt, ta = sum(targets), sum(actual)
    if ta == 0:
        return float("inf")
    return max(abs(a / ta - t / tt) for t, a in zip(targets, actual))


def search_graft(cfg):
    P = sum(cfg["scale_pattern"])
    k_ideal = P / cfg["base_units_per_scale"]
    bx, by = graft_bases(cfg["x"]), graft_bases(cfg["y"])
    lo, hi, steps = k_ideal * 0.5, k_ideal * 2.0, 20000
    seen = {}
    for i in range(steps + 1):
        k = lo + (hi - lo) * i / steps
        cx, cy = graft_counts(cfg["x"], k, P), graft_counts(cfg["y"], k, P)
        key = (tuple(cx), tuple(cy))
        if key in seen:
            continue
        Nx, Ny = sum(cx), sum(cy)
        if not (cfg["min_n"] <= max(Nx, Ny) <= cfg["max_n"]):
            continue
        rx = [c for s, c in zip(cfg["x"], cx) if "base" in s]
        ry = [c for s, c in zip(cfg["y"], cy) if "base" in s]
        if any(v == 0 for v in rx + ry):
            continue
        k_eff = (sum(rx) + sum(ry)) / (sum(bx) + sum(by))
        seen[key] = dict(cx=cx, cy=cy, Nx=Nx, Ny=Ny,
                         err=max(shape_error(bx, rx), shape_error(by, ry)),
                         scale_dev=k_eff / k_ideal - 1)
    tol = cfg.get("scale_size_tolerance", 0.10)
    cands = [c for c in seen.values() if abs(c["scale_dev"]) <= tol]
    cands.sort(key=lambda c: c["err"] + 0.25 * abs(c["scale_dev"])
               + (0.05 if c["Nx"] != c["Ny"] else 0))
    return cands


# --------------------------------------------- proportional-mode search ---

def prop_axis_options(seq, m, cfg):
    """
    All sensible snappings of one axis when every scale is m*P grid units.
    Scale segments are fixed at n*m*P units; each base segment is rounded
    down or up from its ideal size. Proportions are compared on folded
    lengths for 'folded_length' segments. Returns [(counts, err)].
    """
    P = sum(cfg["scale_pattern"])
    L = [seg_length(s) for s in seq]
    tot = sum(L)
    frac = [l / tot for l in L]
    paper = [s["scales"] * m * P if "scales" in s else 0 for s in seq]
    scale_units = sum(f for s, f in zip(seq, folded_counts(seq, paper, cfg))
                      if "scales" in s)
    scale_frac = sum(f for s, f in zip(seq, frac) if "scales" in s)
    N_ideal = scale_units / scale_frac
    base_idx = [i for i, s in enumerate(seq) if "scales" not in s]
    if len(base_idx) > 12:
        choices = [[int(round(frac[i] * N_ideal))] for i in base_idx]
    else:
        choices = [sorted({math.floor(frac[i] * N_ideal),
                           math.ceil(frac[i] * N_ideal)}) for i in base_idx]
    options = []
    for combo in itertools.product(*choices):
        if min(combo, default=1) <= 0:
            continue
        counts = [s["scales"] * m * P if "scales" in s else 0 for s in seq]
        for i, v in zip(base_idx, combo):
            counts[i] = v
        options.append((counts,
                        shape_error(L, folded_counts(seq, counts, cfg))))
    return options


def search_proportional(cfg):
    """
    One family of grids per scale size (m = 1, 2, ...: each scale is m times
    the scale_pattern). Within a family, list the best few roundings of the
    scale-free chunks.
    """
    per_m = cfg.get("options_per_scale_size", 3)
    cands = []
    m = 1
    while True:
        ox = prop_axis_options(cfg["x"], m, cfg)
        oy = prop_axis_options(cfg["y"], m, cfg)
        if not ox or not oy:
            break
        if min(sum(c) for c, _ in ox) > cfg["max_n"]:
            break
        family = []
        for cx, ex in ox:
            for cy, ey in oy:
                N = sum(cx)
                if N != sum(cy) or not (cfg["min_n"] <= N <= cfg["max_n"]):
                    continue
                family.append(dict(cx=cx, cy=cy, Nx=N, Ny=N,
                                   err=max(ex, ey)))
        family.sort(key=lambda c: c["err"])
        cands += family[:per_m]
        m += 1
    return cands


# --------------------------------------------------------- dividing in N ---

def division_reference(N):
    """How to find an N-division reference with diagonals + binary folds."""
    if N & (N - 1) == 0:
        return ("N = %d is a power of two: just keep halving "
                "(%d halvings)." % (N, int(math.log2(N))))
    two_k = 1 << (N.bit_length() - 1)   # largest power of 2 below N
    p = N - two_k
    g = math.gcd(p, two_k)
    return (
        "N = %d = %d + %d. Pinch a mark on the right edge at height %d/%d "
        "(binary folds only). Crease a line from the bottom-left corner "
        "to that mark; where it crosses the diagonal running from top-left to "
        "bottom-right, it is exactly %d/%d of the way across. Crease a "
        "vertical line there. The left part is %d grid units wide (divide it "
        "by halving); then fill the remaining %d units by folding edges "
        "and lines to existing lines."
        % (N, two_k, p, p // g, two_k // g, two_k, N, two_k, p))


# ---------------------------------------------------------------- output ---

def write_csv(path, axes, unit_mm):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["axis", "index", "grid_units", "fraction", "mm_from_start",
                    "mm_from_end", "kind"])
        for axis, lines in axes:
            total = lines[-1][0]
            for i, (p, kind) in enumerate(lines):
                w.writerow([axis, i, "%g" % round(p, 6),
                            "%.5f" % (p / total),
                            "%.2f" % (p * unit_mm),
                            "%.2f" % ((total - p) * unit_mm), kind])


def write_svg(path, xl, yl, xs, ys, title, grid_n=None):
    W = 900.0
    margin = 40.0
    tx, ty = xl[-1][0], yl[-1][0]
    s = (W - 2 * margin) / max(tx, ty)
    Wd, Hd = tx * s + 2 * margin, ty * s + 2 * margin + 30

    def X(u):
        return margin + u * s

    def Y(u):
        return margin + 30 + u * s

    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" '
           'viewBox="0 0 %.1f %.1f">' % (Wd, Hd, Wd, Hd),
           '<rect width="100%" height="100%" fill="white"/>',
           '<text x="%.1f" y="24" font-family="sans-serif" font-size="15">%s'
           '</text>' % (margin, title)]
    # shading: scale strips light, crossings darker
    for (a, b, seg) in xs:
        if "scales" in seg:
            out.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" '
                       'fill="#4a90d9" fill-opacity="0.12"/>'
                       % (X(a), Y(0), (b - a) * s, ty * s))
    for (a, b, seg) in ys:
        if "scales" in seg:
            out.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" '
                       'fill="#4a90d9" fill-opacity="0.12"/>'
                       % (X(0), Y(a), tx * s, (b - a) * s))
    # faint uniform grid underlay (every grid unit) when snapped
    if grid_n:
        for i in range(int(round(tx)) + 1):
            out.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" '
                       'stroke="#ddd" stroke-width="0.4"/>'
                       % (X(i), Y(0), X(i), Y(ty)))
        for i in range(int(round(ty)) + 1):
            out.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" '
                       'stroke="#ddd" stroke-width="0.4"/>'
                       % (X(0), Y(i), X(tx), Y(i)))
    style = {"scale": ('#4a90d9', 0.6), "boundary": ('#d0021b', 1.6),
             "edge": ('#000', 2.0)}
    for p, kind in xl:
        c, wdt = style[kind]
        out.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" '
                   'stroke-width="%.1f"/>' % (X(p), Y(0), X(p), Y(ty), c, wdt))
    for p, kind in yl:
        c, wdt = style[kind]
        out.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" '
                   'stroke-width="%.1f"/>' % (X(0), Y(p), X(tx), Y(p), c, wdt))
    # label chunk boundaries with their grid position
    for p, kind in xl:
        if kind == "boundary":
            out.append('<text x="%.2f" y="%.2f" font-family="sans-serif" '
                       'font-size="10" text-anchor="middle" fill="#d0021b">'
                       '%g</text>' % (X(p), Y(0) - 4, round(p, 3)))
    for p, kind in yl:
        if kind == "boundary":
            out.append('<text x="%.2f" y="%.2f" font-family="sans-serif" '
                       'font-size="10" text-anchor="end" fill="#d0021b">'
                       '%g</text>' % (X(0) - 4, Y(p) + 3, round(p, 3)))
    # chunk names along the top and left
    for (a, b, seg), i in zip(xs, range(len(xs))):
        if "name" in seg:
            out.append('<text x="%.2f" y="%.2f" font-family="sans-serif" '
                       'font-size="11" text-anchor="middle" fill="#555">'
                       '%s</text>' % (X((a + b) / 2), Y(0) - 16, seg["name"]))
    out.append('</svg>')
    with open(path, "w") as f:
        f.write("\n".join(out))


def describe_axis(name, cfg, seq, counts, spans, unit_mm):
    P = sum(cfg["scale_pattern"])
    out = ["%s axis:" % name]
    if cfg["mode"] == "proportional":
        L = [seg_length(s) for s in seq]
        fc = folded_counts(seq, counts, cfg)
        tl, tc = sum(L), sum(fc)
        of = "of folded" if any_folded(cfg) else "of paper"
    for i, ((a, b, seg), c) in enumerate(zip(spans, counts)):
        label = "  %-8s %6g u  %7.1f mm" % (seg_name(seg, i), round(b - a, 4),
                                            (b - a) * unit_mm)
        if cfg["mode"] == "proportional":
            label += "   %.4f %s (target %.4f)" % (fc[i] / tc, of, L[i] / tl)
            if "folded_length" in seg:
                label += "   folds to %g u = %.1f mm" % (
                    round(fc[i], 4), fc[i] * unit_mm)
        else:
            if "base" in seg:
                label += "   (base %g)" % seg["base"]
        if "scales" in seg:
            pitch = c / seg["scales"]
            label += "   %d scales x %g u = %.2f mm each" % (
                seg["scales"], round(pitch, 4), pitch * unit_mm)
        out.append(label)
    return out


def emit(cfg, cx, cy, snap, outdir, header):
    pattern = cfg["scale_pattern"]
    xl, xs = layout(cfg["x"], cx, pattern)
    yl, ys = layout(cfg["y"], cy, pattern)
    tx, ty = xl[-1][0], yl[-1][0]
    pitches = scale_pitches(cfg["x"], cx, pattern) + \
        scale_pitches(cfg["y"], cy, pattern)
    if "scale_width_mm" in cfg:
        unit_mm = cfg["scale_width_mm"] / pitches[0]
        paper = unit_mm * max(tx, ty)
    else:
        paper = cfg["paper_size_mm"]
        unit_mm = paper / max(tx, ty)
    os.makedirs(outdir, exist_ok=True)
    lines = [header, ""]
    if abs(tx - ty) > 1e-9:
        lines.append("WARNING: x total (%g) != y total (%g): this is a "
                     "rectangle, not a square. Adjust the segments so both "
                     "axes match." % (tx, ty))
    if max(pitches) - min(pitches) > 1e-9:
        lines.append("WARNING: scale segments have different scale widths "
                     "(%s grid units) - the scales will not all be the same "
                     "size." % ", ".join("%g" % round(p, 4) for p in pitches))
    if "scale_width_mm" in cfg:
        lines.append("Scale width fixed at %g mm -> paper needed: %.1f mm "
                     "square" % (cfg["scale_width_mm"], paper))
    else:
        lines.append("Paper: %g mm square" % paper)
    lines.append("One grid unit = %.3f mm; one scale = %g grid units = "
                 "%.2f mm (pattern %s)"
                 % (unit_mm, round(pitches[0], 4), pitches[0] * unit_mm,
                    pattern))
    if any_folded(cfg):
        F = cfg["folded_units_per_scale"]
        fp = pitches[0] * F / sum(pattern)
        lines.append("Folded: each scale adds %g grid units = %.2f mm to the "
                     "folded body (folded_units_per_scale = %g)"
                     % (round(fp, 4), fp * unit_mm, F))
    nx = sum(s["scales"] for s in cfg["x"] if "scales" in s)
    ny = sum(s["scales"] for s in cfg["y"] if "scales" in s)
    lines.append("Scales: %d across x, %d across y -> %d scale cells where "
                 "the scale strips cross" % (nx, ny, nx * ny))
    lines.append("")
    lines += describe_axis("x", cfg, cfg["x"], cx, xs, unit_mm)
    lines += describe_axis("y", cfg, cfg["y"], cy, ys, unit_mm)
    if snap and abs(tx - ty) < 1e-9:
        lines += ["", "Dividing the square into %d:" % int(tx),
                  "  " + division_reference(int(tx))]
    summary = "\n".join(lines)
    with open(os.path.join(outdir, "summary.txt"), "w") as f:
        f.write(summary + "\n")
    write_csv(os.path.join(outdir, "lines.csv"), (("x", xl), ("y", yl)),
              unit_mm)
    title = ("%d x %d grid" % (tx, ty)) if snap else "exact (unsnapped) layout"
    title += " - red: chunk boundaries, blue: scale lines"
    write_svg(os.path.join(outdir, "grid.svg"), xl, yl, xs, ys, title,
              grid_n=int(tx) if snap else None)
    print(summary)
    print("\nWrote %s/{summary.txt, lines.csv, grid.svg}" % outdir)


def exact_counts(cfg):
    P = sum(cfg["scale_pattern"])
    if cfg["mode"] == "graft":
        k = P / cfg["base_units_per_scale"]
        return ([s["scales"] * P if "scales" in s else k * s["base"]
                 for s in cfg[a]] for a in ("x", "y"))
    # proportional: scale the layout so the first scale segment's scales are
    # exactly one scale_pattern long, so "grid units" still mean something
    F = cfg.get("folded_units_per_scale", P)
    first = next(s for s in cfg["x"] if "scales" in s)
    first_units = first["scales"] * (F if "folded_length" in first else P)
    k = first_units / seg_length(first)   # target-frame units per length

    def paper_units(s):
        u = seg_length(s) * k
        return u * P / F if "folded_length" in s else u
    return ([paper_units(s) for s in cfg[a]] for a in ("x", "y"))


# ------------------------------------------------------------------ main ---

def print_candidates(cfg, cands, top):
    P = sum(cfg["scale_pattern"])
    names = [seg_name(s, i) for i, s in enumerate(cfg["x"])]
    if cfg["mode"] == "proportional":
        print("Grouped by scale size (grid units per scale); within each, "
              "the best roundings of the scale-free chunks. max err = worst "
              "chunk size error, as a fraction of the %s. Segment sizes "
              "are grid units on the paper.\n"
              % ("folded length" if any_folded(cfg) else "paper"))
        hdr = "  #    N   max err  scale(u)  %s" % "  ".join(names)
        print(hdr)
        for i, c in enumerate(cands[:top], 1):
            pitch = scale_pitches(cfg["x"], c["cx"], cfg["scale_pattern"])[0]
            extra = "" if c["cy"] == c["cx"] else "   y: %s" % c["cy"]
            print("%3d %5d   %5.2f%%  %6g    %s%s" % (
                i, c["Nx"], 100 * c["err"], pitch,
                "  ".join("%*d" % (len(n), v) for n, v in zip(names, c["cx"])),
                extra))
    else:
        print("Graft mode. shape err = worst base-segment error as a "
              "fraction of the base; scale dev = how much bigger(+)/"
              "smaller(-) the base is relative to the scales than "
              "requested.\n")
        print("  #    N   shape err  scale dev   x segments (grid units)")
        for i, c in enumerate(cands[:top], 1):
            n = ("%d" % c["Nx"]) if c["Nx"] == c["Ny"] else \
                ("%dx%d" % (c["Nx"], c["Ny"]))
            print("%3d %5s   %6.2f%%   %+6.1f%%    %s%s" % (
                i, n, 100 * c["err"], 100 * c["scale_dev"], c["cx"],
                "" if c["cy"] == c["cx"] else "  y:%s" % c["cy"]))
    print("\nRun again with --pick # or --n N to write the grid files.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--pick", type=int, help="write files for candidate #")
    ap.add_argument("--n", type=int, help="write files for the grid with N")
    ap.add_argument("--exact", action="store_true",
                    help="no snapping; use the exact requested proportions")
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--out", default="out")
    a = ap.parse_args()
    cfg = load_config(a.config)

    if a.exact:
        cx, cy = exact_counts(cfg)
        emit(cfg, cx, cy, False, a.out, "Exact layout (%s mode)" % cfg["mode"])
        return

    if cfg["mode"] == "proportional":
        cands = search_proportional(cfg)
    else:
        cands = search_graft(cfg)
    if not cands:
        if cfg["mode"] == "proportional":
            sys.exit("No candidates: raise max_n, or check that the scale "
                     "chunks on x and y take up matching shares of the paper "
                     "(the scales are the same size on both axes, so x and y "
                     "must come out to the same total).")
        sys.exit("No candidates: try raising max_n or scale_size_tolerance.")

    if a.n is not None:
        pool = [c for c in cands if c["Nx"] == a.n and c["Ny"] == a.n]
        if not pool:
            sys.exit("No grid with N = %d fits these proportions" % a.n)
        chosen = pool[0]
    elif a.pick is not None:
        if not 1 <= a.pick <= len(cands):
            sys.exit("--pick must be between 1 and %d" % len(cands))
        chosen = cands[a.pick - 1]
    else:
        print_candidates(cfg, cands, a.top)
        return

    n = chosen["Nx"] if chosen["Nx"] == chosen["Ny"] else \
        "%dx%d" % (chosen["Nx"], chosen["Ny"])
    if cfg["mode"] == "proportional":
        header = "N = %s (%s mode), worst chunk error %.2f%% of the %s" % (
            n, cfg["mode"], 100 * chosen["err"],
            "folded length" if any_folded(cfg) else "paper")
    else:
        header = ("N = %s (graft mode), base shape error %.2f%%, scale size "
                  "deviation %+.1f%%" % (n, 100 * chosen["err"],
                                         100 * chosen["scale_dev"]))
    emit(cfg, chosen["cx"], chosen["cy"], True, a.out, header)


if __name__ == "__main__":
    main()
