# Koi scale grid calculator

`koi_grid.py` creates a folding grid for an origami base with a grafted scale
pattern, such as Robert Lang's Koi from *Origami Design Secrets* (Fig. 7.29)
folded with lots of scales (his Scaled Koi, Opus 425, uses a 22 × 22 grid of
scales).

You describe the paper as a sequence of chunks (e.g. head / body / tail) and say
how many scales the body holds. Chunks can also hold other square-grid
tessellations. The script finds a single uniform N × N grid on
which every chunk boundary and every scale crease falls exactly on a grid line,
while keeping your proportions as close as possible. It then writes a diagram of
the grid, a table of every crease position in mm, and instructions for dividing
the square into N.

## Requirements

Python 3.6 or newer. No extra packages.

## Quick start

```sh
# 1. List the candidate grids for a config
python3 koi_grid.py koi_folded.json

# 2. Write the files for the one you like (by row number or by N)
python3 koi_grid.py koi_folded.json --pick 1 --out out_folded
python3 koi_grid.py koi_folded.json --n 132 --out out_132
```

> **The example configs use placeholder values.** Replace the head/body/tail
> ratios and `folded_units_per_scale` with your own before folding anything.

| Example config | Shows |
|---|---|
| `koi_folded.json` | Scales, with the head : body : tail ratio fixed on the folded model |
| `koi_proportional.json` | Scales, with the ratio fixed on the flat paper |
| `koi_tessellation.json` | Scales on the body plus a second tessellation on the tail |
| `example_config.json` | Graft mode |

## Writing a config

A config is a JSON file with some top-level settings and one list of chunks for
each axis of the square:

```json
{
  "paper_size_mm": 600,
  "scale_pattern": [1, 0.5],
  "folded_units_per_scale": 1,
  "x": [
    {"name": "head", "length": 0.25},
    {"name": "body", "folded_length": 0.50, "scales": 20},
    {"name": "tail", "length": 0.25}
  ]
}
```

- **`x`** lists the chunks from left to right; **`y`** lists them from top to
  bottom. Leave `y` out to use a copy of `x`.
- **Top-level settings** (paper size, the default tessellation, grid limits)
  are listed in the [config reference](#config-reference).
- JSON has no comments. Put notes in keys starting with `_`, such as
  `"_comment"`, which are ignored.

### Two kinds of chunk

- A **plain chunk** has no creases of its own, e.g. the head or tail, or a
  piece of the base: `{"name": "head", "length": 0.25}`.
- A **tessellated chunk** is filled with a whole number of repeats of a
  tessellation. Give the count as `"scales"` or `"repeats"` (they mean the
  same thing): `{"name": "body", "folded_length": 0.50, "scales": 20}`.
  Every axis needs at least one tessellated chunk.

### How big each chunk is

| Key | On | Meaning |
|---|---|---|
| `length` | any chunk | Its share of the **flat paper**. Shares are ratios and don't need to add up to 1. |
| `folded_length` | tessellated chunks | Its share **after its pleats are folded flat**. Plain chunks don't fold down, so the `length` of a plain chunk is compared with the `folded_length` of the others. |
| `base` | plain chunks, graft mode | Its length in base units. The tessellated chunks then have no length (see [graft mode](#3-graft-mode)). |

The script picks the mode from the tessellated chunks. If they have a `length`
or `folded_length`, it uses proportional mode; if none do, it uses graft mode.
Either every tessellated chunk has one of the two, or none does. Mixing
`length` and `folded_length` across chunks is fine.

### Which tessellation a chunk uses

A tessellated chunk uses the top-level `scale_pattern` and
`folded_units_per_scale`, unless it has its own `pattern` and `folded_units`:

- **`pattern`** is the spacing between successive creases of one repeat along
  that axis, in grid units. It can include fractions. Its sum is the width of
  one repeat on the paper. `[1, 0.5]` means a crease, 1 unit, a crease, ½ unit,
  then the next repeat.
- **`folded_units`** is how much one repeat adds to the chunk's length once
  folded. It must be more than 0 and no more than the sum of the pattern, and
  it's needed whenever the chunk uses `folded_length`.

Chunks with the same pattern and folded units count as one tessellation and are
always kept the same size. See
[other square-grid tessellations](#other-square-grid-tessellations).

### Rules the script checks

If the config breaks any of these rules, the script stops with a message saying
what's wrong (and which chunk, where it applies):

- counts (`scales` / `repeats`) are whole numbers of at least 1
- a chunk has `scales` or `repeats`, not both, and `length` or
  `folded_length`, not both
- `folded_length`, `pattern` and `folded_units` only go on tessellated chunks
- `paper_size_mm` and `scale_width_mm` aren't both set
- in graft mode, plain chunks use `base`, and each tessellated block
  (count × pattern sum) is a whole number of grid units

## Modes

The mode is picked automatically from how the scale chunks are written.

### 1. Folded proportional mode (recommended)

Fix the head : body : tail ratio **as it looks in the folded model**, where the
body is measured after the scale pleats are collapsed (its length in the base
without scales).

```json
{
  "paper_size_mm": 600,
  "scale_pattern": [1, 1],
  "folded_units_per_scale": 1,
  "x": [
    {"name": "head", "length": 0.30},
    {"name": "body", "folded_length": 0.50, "scales": 22},
    {"name": "tail", "length": 0.20}
  ]
}
```

`folded_units_per_scale` is how many grid units of length one scale adds to the
folded body. Each scale uses `sum(scale_pattern)` units of paper; the rest are
hidden inside the pleats. The best way to find this number is to fold a small
test strip of a few scales and measure it.

**Pleats thinner than a grid unit.** `scale_pattern` can contain fractions.
For pleats ½ unit thick, spaced 1 unit apart:

```json
"scale_pattern": [1, 0.5],
"folded_units_per_scale": 1
```

Each scale then uses 1.5 units of paper and 1 unit of folded length, so the
body folds to ⅔ of its paper width (with `[1, 1]` it would be ½). Chunk
boundaries still land on whole grid lines. The pleat creases that fall between
grid lines are dashed in `grid.svg`, marked `scale-sub` in `lines.csv`, and the
summary says how finely to divide those grid cells. If a scale count would
leave the body a fractional number of units long (e.g. 5 × 1.5 = 7.5), that
grid is skipped and the next scale size is used instead.

See `koi_folded.json`.

### 2. Paper proportional mode

The same as above, but use `"length"` on the scale chunk, so the ratio is
measured on the **flat paper** instead of the folded model.

```json
{"name": "body", "length": 0.50, "scales": 22}
```

See `koi_proportional.json`.

### 3. Graft mode

The original approach: cut the scale-free base at some points and insert blocks
of scales there. Base pieces use `"base"` (in any unit you like), scale blocks
have no length, and `base_units_per_scale` sets how big one scale is compared to
the base. The body's length then depends on the scale size and count.

```json
{
  "base_units_per_scale": 0.05,
  "x": [ {"base": 0.30}, {"scales": 22}, {"base": 0.70} ]
}
```

Other tessellations can be grafted in the same way, e.g.
`{"repeats": 3, "pattern": [2, 1, 1]}`. Each block is its count × its pattern
sum grid units long. `base_units_per_scale` refers to the first tessellated
chunk.

See `example_config.json`.

## Other square-grid tessellations

Scales are just one tessellation. Any chunk can hold a different square-grid
tessellation by giving it its own `pattern` (the crease offsets of one repeat
along that axis) and `folded_units` (how much one repeat adds to the folded
length). These override the global `scale_pattern` and
`folded_units_per_scale`. You can write `repeats` instead of `scales`
(see `koi_tessellation.json`):

```json
"x": [
  {"name": "head", "length": 0.25},
  {"name": "body", "folded_length": 0.50, "scales": 20},
  {"name": "tail", "folded_length": 0.25, "repeats": 3,
   "pattern": [2, 1, 1], "folded_units": 2}
]
```

- **Same pattern, same size.** Chunks with the same `pattern` and
  `folded_units` are one tessellation: their repeats are always identical,
  on both axes (e.g. a body split around a base line).
- **Different patterns are sized independently.** The first tessellation
  (the first tessellated chunk in `x`) sets the grid. The `scale(u)` and mm
  columns in the candidate list, `scale_width_mm` and `base_units_per_scale`
  all refer to it. Every other tessellation is scaled to fit its own chunk, in steps
  of ½ its pattern by default, so some of its creases may fall on half grid
  lines (dashed, like ½-unit pleats). Set `tessellation_subdivision` to 1 to
  allow only whole steps (coarser fit, every crease on the grid), or higher
  for a closer fit with finer creases.
- **Different patterns on x and y.** Give `x` and `y` their own chunk lists
  with different patterns, for tessellations that aren't symmetric.

The script only places the horizontal and vertical lines. Diagonals and the
mountain/valley assignment come from the tessellation's own crease pattern,
using these lines as references. Tessellations whose vertices aren't on a
square grid (e.g. some 22.5° twists) can't be snapped; use `--exact` for
measured positions instead.

## Config reference

| Key | Default | Meaning |
|---|---|---|
| `x` | *(required)* | List of chunks along the x axis, in order (see below). |
| `y` | copy of `x` | List of chunks along the y axis. If the fish lies along the paper's diagonal, use the same list for both axes. |
| `paper_size_mm` | 500 | Side of the square paper. The scale width is reported. |
| `scale_width_mm` | – | Use **instead of** `paper_size_mm` to fix how wide one repeat of the first tessellation is on the paper (its full pattern, before folding). The paper size needed is reported. |
| `scale_pattern` | `[1, 1]` | Default `pattern` for tessellated chunks that don't have their own. `[1, 1]` = two equal parts; `[1, 2]` = a 1-unit part and a 2-unit part; `[1, 0.5]` = a 1-unit part and a ½-unit pleat (see above). |
| `folded_units_per_scale` | – | Default `folded_units` for tessellated chunks that don't have their own. Needed if a `folded_length` chunk has no `folded_units`. |
| `max_n` / `min_n` | 160 / 8 | Largest / smallest grid to consider. |
| `options_per_scale_size` | 3 | Proportional modes: how many roundings to list for each size of the first tessellation. |
| `base_units_per_scale` | – | Graft mode only: length of one repeat of the first tessellation, in base units. |
| `scale_size_tolerance` | 0.10 | Graft mode only: how far (±10%) the scale size may drift from the requested one. |
| `tessellation_subdivision` | 2 | When there are several tessellations: every one after the first is sized in steps of 1/this of its pattern. |

Keys starting with `_` (e.g. `_comment`) are ignored, so you can use them for notes.

### Chunk keys

| Key | On | Meaning |
|---|---|---|
| `name` | any chunk | Label used in the output. |
| `length` | any chunk | Share of the flat paper. |
| `folded_length` | tessellated | Share after folding. |
| `base` | plain, graft mode | Length in base units. |
| `scales` / `repeats` | tessellated | Number of repeats (whole number, at least 1). Use one or the other. |
| `pattern` | tessellated | This chunk's own crease spacing; overrides `scale_pattern`. |
| `folded_units` | tessellated | This chunk's own folded units per repeat; overrides `folded_units_per_scale`. |

See [Writing a config](#writing-a-config) for how these fit together. You can split chunks
further, e.g. to keep a line from the base inside the head, or to split the body
around a line: `{"scales": 11, ...}, {"length": 0.05}, {"scales": 11, ...}`.
Chunks with the same pattern are always kept the same size.

## Reading the candidate list

```
  #    N   max err  scale(u)  scale mm  folded mm  head  body  tail
  1    66    0.45%       2      18.18       9.09     13    44     9
  4   132    0.45%       4      18.18       9.09     26    88    18
```

- **N**: the grid is N × N.
- **max err**: the worst difference between a chunk's actual and target share,
  as a fraction of the (folded or flat) length. 0.45% of a 600 mm sheet is
  under 3 mm.
- **scale(u)**: grid units per scale (per repeat of the first tessellation)
  on the paper. Every scale is identical and
  the body is always a whole number of scales, so only the scale-free chunks
  get rounded. Bigger scale sizes mean finer grids and more accurate ratios,
  but many more creases.
- **scale mm**: length of one scale on the flat paper (its full
  `scale_pattern`). **folded mm** (only with `folded_length`): how much one
  scale adds to the folded body. If you set `scale_width_mm` instead of
  `paper_size_mm`, this column shows the **paper mm** needed.
- The remaining columns are each chunk's size in grid units on the paper.

## Output files

Written to the `--out` directory (default `out/`):

| File | Contents |
|---|---|
| `grid.svg` | Diagram of the square: black = paper edges, red = chunk boundaries (labelled with their grid position), blue = tessellation creases (dashed if between grid lines), shaded = tessellation strips (darker where they cross; a different colour for each tessellation). |
| `lines.csv` | Every line on each axis: grid index, fraction of the paper, mm from each edge, and its kind: `edge`, `boundary`, `scale` (a crease of any tessellation), or `scale-sub` (one that falls between grid lines). |
| `summary.txt` | Paper size, the size of one repeat of each tessellation (flat and folded), each chunk's size and share compared with the target, and how to divide the square into N. |

## Dividing the square into N

If N is a power of two, keep halving. Otherwise the summary gives a single
reference fold. For example, for N = 88 = 64 + 24:

1. Pinch a mark on the right edge at 24/64 = 3/8 of the way up (binary folds only).
2. Crease a line from the bottom-left corner to that mark.
3. Where it crosses the diagonal running from top-left to bottom-right is
   exactly 64/88 of the way across. Crease a vertical line there.
4. Halve the left part down to single units, then fill in the remaining 24
   units by folding edges and lines to existing lines.

## Command-line options

| Option | Meaning |
|---|---|
| `--pick K` | Write the files for row K of the candidate list. |
| `--n N` | Write the files for the best listed grid with this N. |
| `--exact` | Skip grid snapping and use your exact proportions. Useful if you'd rather measure creases in mm than fold a grid. |
| `--top K` | How many candidates to list (default 15). |
| `--out DIR` | Where to write the files (default `out`). |

## Caveats

- The script only lays out the grid lines. It doesn't assign mountain and valley
  folds or draw the rest of the crease pattern; take those from the book.
- None of the example values are Lang's. Measure your base and scale pattern
  and put those numbers in your config.
- Only square-grid tessellations are supported. Triangular grids would need
  three directions of lines, which this layout doesn't handle.
