#!/usr/bin/env python3
"""
koi_grid.py - compute a folding grid for a base with grafted tessellations
(e.g. scales).

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
  of the first scale segment should be, in base units).

In both modes each scale occupies `scale_pattern` grid units (default [1, 1]:
e.g. one pleat unit + one visible unit), possibly multiplied by a whole number
so the proportions can be matched on a finer grid.

OTHER SQUARE-GRID TESSELLATIONS: any scale segment can carry its own
"pattern" (crease offsets of one repeat of the tessellation along this axis)
and "folded_units" (how much one repeat adds to the folded length), which
override scale_pattern / folded_units_per_scale. "repeats" may be used
instead of "scales". Segments with the same pattern and folded_units are
kept identical in size; different tessellations are each sized to fit their
own chunk.
  {"name": "tail", "length": 0.2, "repeats": 4, "pattern": [2, 1, 1]}

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
import copy
import csv
import itertools
import json
import math
import os
import sys


# ---------------------------------------------------------------- config ---
#
# load_config() annotates every tessellated segment (one with "scales") with:
#   _pattern  crease offsets of one repeat, _P = sum(_pattern)
#   _F        folded units per repeat (None if not given)
#   _group    (pattern, folded units): segments in a group stay identical
#   _word     "scales" or "repeats", for the output

def seg_length(seg):
    """Target proportional length of a segment (proportional mode)."""
    for key in ("folded_length", "length", "base"):
        if key in seg:
            return seg[key]


def fold(seg, c):
    """Length of a segment of c paper units once its pleats are collapsed.
    Only 'folded_length' segments are measured folded."""
    return c * seg["_F"] / seg["_P"] if "folded_length" in seg else c


def unfold(seg, f):
    """Inverse of fold()."""
    return f * seg["_P"] / seg["_F"] if "folded_length" in seg else f


def folded_counts(seq, counts):
    return [fold(s, c) for s, c in zip(seq, counts)]


def is_whole(v, tol=1e-9):
    return abs(v - round(v)) < tol


def any_folded(cfg):
    return any("folded_length" in s for a in ("x", "y") for s in cfg[a])


def tess_segs(cfg, axes=("x", "y")):
    return [s for a in axes for s in cfg[a] if "scales" in s]


def groups(cfg):
    """Distinct tessellations, in order of first appearance (x, then y)."""
    out = []
    for s in tess_segs(cfg):
        if s["_group"] not in out:
            out.append(s["_group"])
    return out


def group_label(cfg, g):
    names, axes = [], []
    for axis in ("x", "y"):
        for i, s in enumerate(cfg[axis]):
            if "scales" in s and s["_group"] == g:
                if seg_name(s, i) not in names:
                    names.append(seg_name(s, i))
                if axis not in axes:
                    axes.append(axis)
    label = "/".join(names)
    return label if len(axes) == 2 else "%s (%s only)" % (label, axes[0])


def seg_name(seg, i):
    if "name" in seg:
        return seg["name"]
    return ("scales%d" % i) if "scales" in seg else ("base%d" % i)


def show(seg):
    """A segment as the user wrote it, for error messages."""
    return {k: v for k, v in seg.items() if not k.startswith("_")}


def load_config(path):
    with open(path) as f:
        cfg = json.load(f)
    cfg.setdefault("scale_pattern", [1, 1])
    cfg.setdefault("max_n", 160)
    cfg.setdefault("min_n", 8)
    if "y" not in cfg:
        cfg["y"] = copy.deepcopy(cfg["x"])
    if "paper_size_mm" in cfg and "scale_width_mm" in cfg:
        sys.exit("set either 'paper_size_mm' or 'scale_width_mm', not both")
    if "paper_size_mm" not in cfg and "scale_width_mm" not in cfg:
        cfg["paper_size_mm"] = 500

    scale_segs = []
    for axis in ("x", "y"):
        for seg in cfg[axis]:
            if "repeats" in seg:
                if "scales" in seg:
                    sys.exit("give 'scales' or 'repeats', not both: %r"
                             % show(seg))
                seg["scales"] = seg["repeats"]
                seg["_word"] = "repeats"
            if "scales" in seg:
                seg.setdefault("_word", "scales")
                if "base" in seg:
                    sys.exit("segment has both 'base' and 'scales': %r"
                             % show(seg))
                if seg["scales"] < 1 or not is_whole(seg["scales"]):
                    sys.exit("scale/repeat count must be a whole number "
                             ">= 1: %r" % show(seg))
                seg["scales"] = int(round(seg["scales"]))
                pattern = seg.get("pattern", cfg["scale_pattern"])
                if not pattern or any(not isinstance(u, (int, float))
                                      or u <= 0 for u in pattern):
                    sys.exit("a pattern must be a list of positive numbers: "
                             "%r" % show(seg))
                seg["_pattern"] = list(pattern)
                seg["_P"] = sum(pattern)
                seg["_F"] = seg.get("folded_units",
                                    cfg.get("folded_units_per_scale"))
                seg["_group"] = (tuple(pattern), seg["_F"])
                scale_segs.append(seg)
            else:
                for key in ("folded_length", "pattern", "folded_units"):
                    if key in seg:
                        sys.exit("'%s' only applies to scale/repeat segments: "
                                 "%r" % (key, show(seg)))
                if "base" not in seg and "length" not in seg:
                    sys.exit("segment needs 'base', 'length' or 'scales': %r"
                             % show(seg))
        if not any("scales" in s for s in cfg[axis]):
            sys.exit("axis '%s' has no scale segment" % axis)
    for s in scale_segs:
        if "length" in s and "folded_length" in s:
            sys.exit("give a scale segment 'length' or 'folded_length', "
                     "not both: %r" % show(s))
        if "folded_length" in s:
            if s["_F"] is None:
                sys.exit("'folded_length' needs 'folded_units_per_scale' (or "
                         "'folded_units' on the segment): how many grid units "
                         "of length one repeat adds to the folded chunk (out "
                         "of the %g units it uses on the paper): %r"
                         % (s["_P"], show(s)))
            if not 0 < s["_F"] <= s["_P"]:
                sys.exit("folded units must be > 0 and <= %g (sum of the "
                         "pattern): %r" % (s["_P"], show(s)))
    with_len = [s for s in scale_segs
                if "length" in s or "folded_length" in s]
    if with_len and len(with_len) != len(scale_segs):
        sys.exit("either every scale segment has a 'length'/'folded_length' "
                 "(proportional mode) or none does (graft mode)")
    cfg["mode"] = "proportional" if with_len else "graft"
    if cfg["mode"] == "graft":
        if "base_units_per_scale" not in cfg:
            sys.exit("graft mode needs 'base_units_per_scale': how long one "
                     "scale should be, in base units")
        for axis in ("x", "y"):
            for seg in cfg[axis]:
                if "length" in seg and "scales" not in seg:
                    sys.exit("graft mode: use 'base' (not 'length') for "
                             "base segments: %r" % show(seg))
                if "scales" in seg and not is_whole(seg["scales"] * seg["_P"]):
                    sys.exit("graft mode: %d scales x %g units = %g is not a "
                             "whole number of grid units; change the scale "
                             "count or pattern"
                             % (seg["scales"], seg["_P"],
                                seg["scales"] * seg["_P"]))
    return cfg


# ---------------------------------------------------------------- layout ---

def layout(seq, counts):
    """
    Place every line on one axis. counts[i] = length of segment i in grid
    units (ints when snapped, floats for exact layouts). Scale segments are
    divided evenly into repeats, each split according to its pattern.
    Returns (lines, spans): lines = [(pos, kind)] with kind in
    edge/boundary/scale, spans = [(start, end, seg)].
    """
    pos = 0.0
    lines = [(0.0, "edge")]
    spans = []
    for seg, c in zip(seq, counts):
        start = pos
        if "scales" in seg:
            unit = c / (seg["scales"] * seg["_P"])
            p = start
            for _ in range(seg["scales"]):
                for u in seg["_pattern"]:
                    p += u * unit
                    lines.append((p, "scale"))
            lines.pop()             # the last one is the segment boundary
        pos = start + c
        lines.append((pos, "boundary"))
        spans.append((start, pos, seg))
    lines[-1] = (pos, "edge")
    return lines, spans


def scale_pitches(seq, counts):
    """Grid units per repeat for every scale segment on an axis."""
    return [c / s["scales"] for s, c in zip(seq, counts) if "scales" in s]


# ---------------------------------------------------- graft-mode search ---

def graft_counts(seq, k):
    return [s["scales"] * s["_P"] if "scales" in s
            else int(round(k * s["base"])) for s in seq]


def graft_bases(seq):
    return [s["base"] for s in seq if "base" in s]


def shape_error(targets, actual):
    """Largest difference between target and actual fractions."""
    tt, ta = sum(targets), sum(actual)
    if ta == 0:
        return float("inf")
    return max(abs(a / ta - t / tt) for t, a in zip(targets, actual))


def search_graft(cfg):
    P = tess_segs(cfg)[0]["_P"]
    k_ideal = P / cfg["base_units_per_scale"]
    bx, by = graft_bases(cfg["x"]), graft_bases(cfg["y"])
    lo, hi, steps = k_ideal * 0.5, k_ideal * 2.0, 20000
    seen = {}
    for i in range(steps + 1):
        k = lo + (hi - lo) * i / steps
        cx, cy = graft_counts(cfg["x"], k), graft_counts(cfg["y"], k)
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

def tess_units(seg, m):
    """Paper grid units of a tessellated segment whose repeats are m times
    its pattern."""
    return seg["scales"] * m * seg["_P"]


def ideal_total(seq, primary, m):
    """Ideal grid size (in the target frame) implied by the primary
    tessellation at multiplier m, or None if it isn't on this axis."""
    L = [seg_length(s) for s in seq]
    tot = sum(L)
    idx = [i for i, s in enumerate(seq)
           if "scales" in s and s["_group"] == primary]
    if not idx:
        return None
    units = sum(fold(seq[i], tess_units(seq[i], m)) for i in idx)
    return units / (sum(L[i] for i in idx) / tot)


def prop_axis_options(seq, primary, m, N_ideal, subdiv):
    """
    All sensible snappings of one axis. Primary-tessellation segments use
    multiplier m; every other tessellation gets the multipliers nearest its
    ideal size, in steps of 1/subdiv (so its creases may fall on 1/subdiv
    grid lines), and each scale-free segment is rounded down or up. Proportions are compared on folded lengths for 'folded_length'
    segments. A tessellated segment must be a whole number of grid units so
    its ends land on the grid. Returns [(counts, err, mults)], or None if
    the primary tessellation can't use multiplier m.
    """
    L = [seg_length(s) for s in seq]
    tot = sum(L)
    frac = [l / tot for l in L]
    for s in seq:
        if "scales" in s and s["_group"] == primary \
                and not is_whole(tess_units(s, m)):
            return None
    others = []
    for s in seq:
        if "scales" in s and s["_group"] != primary \
                and s["_group"] not in others:
            others.append(s["_group"])
    mult_choices = []
    for g in others:
        idx = [i for i, s in enumerate(seq)
               if "scales" in s and s["_group"] == g]
        target = sum(unfold(seq[i], frac[i] * N_ideal) for i in idx)
        ideal = target / sum(tess_units(seq[i], 1) for i in idx)
        ks = [j / subdiv for j in range(max(1, math.floor(ideal * subdiv) - 2),
                                        math.ceil(ideal * subdiv) + 3)]
        ks = [k for k in ks
              if all(is_whole(tess_units(seq[i], k)) for i in idx)]
        ks.sort(key=lambda k: abs(k - ideal))
        if not ks:
            return []
        mult_choices.append(ks[:2])
    base_idx = [i for i, s in enumerate(seq) if "scales" not in s]
    if len(base_idx) > 12:
        choices = [[int(round(frac[i] * N_ideal))] for i in base_idx]
    else:
        choices = [sorted({math.floor(frac[i] * N_ideal),
                           math.ceil(frac[i] * N_ideal)}) for i in base_idx]
    options = []
    for mcombo in itertools.product(*mult_choices):
        mults = dict(zip(others, mcombo))
        mults[primary] = m
        for combo in itertools.product(*choices):
            if min(combo, default=1) <= 0:
                continue
            counts = [int(round(tess_units(s, mults[s["_group"]])))
                      if "scales" in s else 0 for s in seq]
            for i, v in zip(base_idx, combo):
                counts[i] = v
            options.append((counts,
                            shape_error(L, folded_counts(seq, counts)),
                            mults))
    return options


def search_proportional(cfg):
    """
    One family of grids per size of the first tessellation (m = 1, 2, ...:
    each repeat is m times its pattern). Within a family, list the best few
    roundings of everything else. A tessellation that appears on both axes
    uses the same multiplier on both, so its cells stay the same size.
    """
    per_m = cfg.get("options_per_scale_size", 3)
    subdiv = cfg.get("tessellation_subdivision", 2)
    primary = groups(cfg)[0]
    min_units = min(sum(tess_units(s, 1) for s in cfg[a]
                        if "scales" in s and s["_group"] == primary)
                    or float("inf") for a in ("x", "y"))
    cands = []
    m = 0
    while True:
        m += 1
        if m * min_units > cfg["max_n"]:
            break
        nx = ideal_total(cfg["x"], primary, m)
        ny = ideal_total(cfg["y"], primary, m)
        ox = prop_axis_options(cfg["x"], primary, m, nx, subdiv)
        oy = prop_axis_options(cfg["y"], primary, m, ny if ny else nx, subdiv)
        if ox is None or oy is None:
            continue
        if not ox or not oy:
            continue
        if min(sum(c) for c, _, _ in ox) > cfg["max_n"]:
            break
        family = []
        for cx, ex, mx in ox:
            for cy, ey, my in oy:
                N = sum(cx)
                if N != sum(cy) or not (cfg["min_n"] <= N <= cfg["max_n"]):
                    continue
                if any(mx[g] != my[g] for g in mx if g in my):
                    continue
                family.append(dict(cx=cx, cy=cy, Nx=N, Ny=N,
                                   err=max(ex, ey)))
        family.sort(key=lambda c: c["err"])
        cands += family[:per_m]
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

GROUP_COLOURS = ["#4a90d9", "#e8a33d", "#50b36b", "#9b59b6", "#d35d8a"]


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


def write_svg(path, xl, yl, xs, ys, title, colour, grid_n=None):
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
    # shading: tessellated strips light, crossings darker; one colour per
    # tessellation
    for (a, b, seg) in xs:
        if "scales" in seg:
            out.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" '
                       'fill="%s" fill-opacity="0.12"/>'
                       % (X(a), Y(0), (b - a) * s, ty * s, colour(seg)))
    for (a, b, seg) in ys:
        if "scales" in seg:
            out.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" '
                       'fill="%s" fill-opacity="0.12"/>'
                       % (X(0), Y(a), tx * s, (b - a) * s, colour(seg)))
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
    style = {"scale": ('#4a90d9', 0.6, ''),
             "scale-sub": ('#4a90d9', 0.6, ' stroke-dasharray="3 2"'),
             "boundary": ('#d0021b', 1.6, ''), "edge": ('#000', 2.0, '')}
    for p, kind in xl:
        c, wdt, dash = style[kind]
        out.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" '
                   'stroke-width="%.1f"%s/>'
                   % (X(p), Y(0), X(p), Y(ty), c, wdt, dash))
    for p, kind in yl:
        c, wdt, dash = style[kind]
        out.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" '
                   'stroke-width="%.1f"%s/>'
                   % (X(0), Y(p), X(tx), Y(p), c, wdt, dash))
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
    multi = len(groups(cfg)) > 1
    out = ["%s axis:" % name]
    if cfg["mode"] == "proportional":
        L = [seg_length(s) for s in seq]
        fc = folded_counts(seq, counts)
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
            label += "   %d %s x %g u = %.2f mm each" % (
                seg["scales"], seg["_word"], round(pitch, 4), pitch * unit_mm)
            if multi:
                label += " (pattern %s)" % seg["_pattern"]
        out.append(label)
    return out


def emit(cfg, cx, cy, snap, outdir, header):
    xl, xs = layout(cfg["x"], cx)
    yl, ys = layout(cfg["y"], cy)
    if snap:
        xl, yl = mark_off_grid(xl), mark_off_grid(yl)
    tx, ty = xl[-1][0], yl[-1][0]
    # repeat width of every tessellated segment, by tessellation
    by_group = {}
    for seq, counts in ((cfg["x"], cx), (cfg["y"], cy)):
        for s, c in zip(seq, counts):
            if "scales" in s:
                by_group.setdefault(s["_group"], []).append(
                    (s, c / s["scales"]))
    gs = groups(cfg)
    first_seg, first_pitch = by_group[gs[0]][0]
    if "scale_width_mm" in cfg:
        unit_mm = cfg["scale_width_mm"] / first_pitch
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
    for g in gs:
        pitches = [p for _, p in by_group[g]]
        if max(pitches) - min(pitches) > 1e-9:
            lines.append("WARNING: %s segments have different widths "
                         "(%s grid units) - they will not all be the same "
                         "size." % (group_label(cfg, g) if len(gs) > 1
                                    else "scale",
                                    ", ".join("%g" % round(p, 4)
                                              for p in pitches)))
    if "scale_width_mm" in cfg:
        lines.append("Scale width fixed at %g mm -> paper needed: %.1f mm "
                     "square" % (cfg["scale_width_mm"], paper))
    else:
        lines.append("Paper: %g mm square" % paper)
    if len(gs) == 1:
        lines.append("One grid unit = %.3f mm; one scale = %g grid units = "
                     "%.2f mm (pattern %s)"
                     % (unit_mm, round(first_pitch, 4), first_pitch * unit_mm,
                        first_seg["_pattern"]))
        if any_folded(cfg):
            F = first_seg["_F"]
            fp = first_pitch * F / first_seg["_P"]
            lines.append("Folded: each scale adds %g grid units = %.2f mm to "
                         "the folded body (folded_units_per_scale = %g)"
                         % (round(fp, 4), fp * unit_mm, F))
        nx = sum(s["scales"] for s in cfg["x"] if "scales" in s)
        ny = sum(s["scales"] for s in cfg["y"] if "scales" in s)
        lines.append("Scales: %d across x, %d across y -> %d scale cells "
                     "where the scale strips cross" % (nx, ny, nx * ny))
    else:
        lines.append("One grid unit = %.3f mm" % unit_mm)
        lines.append("Tessellations:")
        for g in gs:
            s, pitch = by_group[g][0]
            text = ("  %s: pattern %s, one repeat = %g u = %.2f mm"
                    % (group_label(cfg, g), s["_pattern"], round(pitch, 4),
                       pitch * unit_mm))
            if s["_F"] is not None:
                fp = pitch * s["_F"] / s["_P"]
                text += ", folds to %g u = %.2f mm" % (round(fp, 4),
                                                        fp * unit_mm)
            nx = sum(t["scales"] for t in cfg["x"]
                     if "scales" in t and t["_group"] == g)
            ny = sum(t["scales"] for t in cfg["y"]
                     if "scales" in t and t["_group"] == g)
            text += "; " + ", ".join("%d %s across %s" % (n, s["_word"], a)
                                     for n, a in ((nx, "x"), (ny, "y")) if n)
            lines.append(text)
    off = [p for p, kind in xl + yl if kind == "scale-sub"]
    if off:
        d = next((d for d in range(2, 65)
                  if all(is_whole(p * d, 1e-6) for p in off)), None)
        step = ("1/%d of a grid unit" % d) if d else "fractions of a grid unit"
        lines.append("%d scale creases fall between grid lines, at multiples "
                     "of %s (dashed in grid.svg, 'scale-sub' in lines.csv). "
                     "Fold the %dx%d grid first, then add these by dividing "
                     "the grid cells in the %s strips%s."
                     % (len(off), step, round(tx), round(ty),
                        "tessellation" if len(gs) > 1 else "scale",
                        (" into %d" % d) if d else ""))
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

    def colour(seg):
        return GROUP_COLOURS[gs.index(seg["_group"]) % len(GROUP_COLOURS)]
    write_svg(os.path.join(outdir, "grid.svg"), xl, yl, xs, ys, title,
              colour, grid_n=int(tx) if snap else None)
    print(summary)
    print("\nWrote %s/{summary.txt, lines.csv, grid.svg}" % outdir)


def mark_off_grid(lines):
    """Relabel scale creases that don't sit on a whole grid line."""
    return [(p, "scale-sub" if kind == "scale" and not is_whole(p, 1e-6)
             else kind) for p, kind in lines]


def exact_counts(cfg):
    first = tess_segs(cfg, ("x",))[0]
    if cfg["mode"] == "graft":
        k = first["_P"] / cfg["base_units_per_scale"]
        return ([s["scales"] * s["_P"] if "scales" in s else k * s["base"]
                 for s in cfg[a]] for a in ("x", "y"))
    # proportional: scale the layout so the first scale segment's repeats
    # are exactly one pattern long, so "grid units" still mean something
    first_units = fold(first, first["scales"] * first["_P"])
    k = first_units / seg_length(first)   # target-frame units per length
    return ([unfold(s, seg_length(s) * k) for s in cfg[a]]
            for a in ("x", "y"))


# ------------------------------------------------------------------ main ---

def mm_columns(cfg):
    """Header and per-candidate formatter for the size of one repeat of the
    first tessellation in mm. With scale_width_mm that size is fixed, so the
    paper size is shown."""
    first = tess_segs(cfg)[0]
    folded = any("folded_length" in s for s in tess_segs(cfg)
                 if s["_group"] == first["_group"])
    if "scale_width_mm" in cfg:
        def cols(pitch, N):
            return "%8.1f" % (cfg["scale_width_mm"] / pitch * N)
        return "paper mm", cols
    paper = cfg["paper_size_mm"]

    def cols(pitch, N):
        s = "%8.2f" % (pitch * paper / N)
        if folded:
            s += "  %9.2f" % (pitch * paper / N * first["_F"] / first["_P"])
        return s
    return ("scale mm  folded mm" if folded else "scale mm"), cols


def print_candidates(cfg, cands, top):
    first = tess_segs(cfg)[0]
    mm_hdr, mm_cols = mm_columns(cfg)
    names = [seg_name(s, i) for i, s in enumerate(cfg["x"])]
    gs = groups(cfg)
    if len(gs) > 1:
        print("%d tessellations. The scale(u) and mm columns are for the "
              "first one (%s).\n" % (len(gs), group_label(cfg, gs[0])))
    if cfg["mode"] == "proportional":
        print("Grouped by scale size (grid units per scale); within each, "
              "the best roundings of the scale-free chunks. max err = worst "
              "chunk size error, as a fraction of the %s. Segment sizes "
              "are grid units on the paper.\n"
              % ("folded length" if any_folded(cfg) else "paper"))
        hdr = "  #    N   max err  scale(u)  %s  %s" % (mm_hdr,
                                                     "  ".join(names))
        print(hdr)
        for i, c in enumerate(cands[:top], 1):
            pitch = next(cnt / s["scales"] for s, cnt in zip(cfg["x"], c["cx"])
                         if "scales" in s and s["_group"] == first["_group"])
            extra = "" if c["cy"] == c["cx"] else "   y: %s" % c["cy"]
            print("%3d %5d   %5.2f%%  %6g  %s    %s%s" % (
                i, c["Nx"], 100 * c["err"], pitch,
                mm_cols(pitch, max(c["Nx"], c["Ny"])),
                "  ".join("%*d" % (len(n), v) for n, v in zip(names, c["cx"])),
                extra))
    else:
        print("Graft mode. shape err = worst base-segment error as a "
              "fraction of the base; scale dev = how much bigger(+)/"
              "smaller(-) the base is relative to the scales than "
              "requested.\n")
        print("  #    N   shape err  scale dev  %s   x segments (grid units)"
              % mm_hdr)
        for i, c in enumerate(cands[:top], 1):
            n = ("%d" % c["Nx"]) if c["Nx"] == c["Ny"] else \
                ("%dx%d" % (c["Nx"], c["Ny"]))
            print("%3d %5s   %6.2f%%   %+6.1f%%  %s    %s%s" % (
                i, n, 100 * c["err"], 100 * c["scale_dev"],
                mm_cols(first["_P"], max(c["Nx"], c["Ny"])), c["cx"],
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
