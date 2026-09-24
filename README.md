# Koi scale grid calculator

`koi_grid.py` works out a folding grid for an origami base with a grafted scale
pattern, such as Robert Lang's Koi from *Origami Design Secrets* (Fig. 7.29)
folded with lots of scales (his Scaled Koi, Opus 425, uses a 22 × 22 grid of
scales).

You describe the paper as a sequence of chunks (e.g. head / body / tail) and say
how many scales the body holds. The script finds a single uniform N × N grid on
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

See `example_config.json`.

## Config reference

| Key | Default | Meaning |
|---|---|---|
| `x` | *(required)* | List of chunks along the x axis, in order (see below). |
| `y` | same as `x` | List of chunks along the y axis. If the fish lies along the paper's diagonal, use the same list for both axes. |
| `paper_size_mm` | 500 | Side of the square paper. The scale width is reported. |
| `scale_width_mm` | – | Use **instead of** `paper_size_mm` to fix how wide one scale is on the paper (its full `scale_pattern`, before folding). The paper size needed is reported. |
| `scale_pattern` | `[1, 1]` | Grid units per scale, split into the creases within one scale. `[1, 1]` = two equal parts; `[1, 2]` = a 1-unit part and a 2-unit part; `[1, 0.5]` = a 1-unit part and a ½-unit pleat (see above). |
| `folded_units_per_scale` | – | Required when any chunk uses `folded_length`. See above. |
| `max_n` / `min_n` | 160 / 8 | Largest / smallest grid to consider. |
| `options_per_scale_size` | 3 | Proportional modes: how many roundings to list for each scale size. |
| `base_units_per_scale` | – | Graft mode only: length of one scale in base units. |
| `scale_size_tolerance` | 0.10 | Graft mode only: how far (±10%) the scale size may drift from the requested one. |

Keys starting with `_` (e.g. `_comment`) are ignored, so you can use them for notes.

### Chunks

| Chunk | Meaning |
|---|---|
| `{"length": L}` | A scale-free chunk, `L` = its share (lengths are ratios and need not add up to 1). |
| `{"length": L, "scales": n}` | `n` scales filling a share `L` of the paper. |
| `{"folded_length": L, "scales": n}` | `n` scales whose **folded** length is share `L`. |
| `{"base": L}` / `{"scales": n}` | Graft mode pieces. |

Any chunk can have a `"name"`, which is used in the output. You can split chunks
further, e.g. to keep a line from the base inside the head, or to split the body
around a line: `{"scales": 11, ...}, {"length": 0.05}, {"scales": 11, ...}`.
Keep the scales on every chunk the same size, or the output will warn you.

## Reading the candidate list

```
  #    N   max err  scale(u)  head  body  tail
  1    66    0.45%       2      13    44     9
  4   132    0.45%       4      26    88    18
```

- **N**: the grid is N × N.
- **max err**: the worst difference between a chunk's actual and target share,
  as a fraction of the (folded or flat) length. 0.45% of a 600 mm sheet is
  under 3 mm.
- **scale(u)**: grid units per scale on the paper. Every scale is identical and
  the body is always a whole number of scales, so only the scale-free chunks
  get rounded. Bigger scale sizes mean finer grids and more accurate ratios,
  but many more creases.
- The remaining columns are each chunk's size in grid units on the paper.

## Output files

Written to the `--out` directory (default `out/`):

| File | Contents |
|---|---|
| `grid.svg` | Diagram of the square: black = paper edges, red = chunk boundaries (labelled with their grid position), blue = scale creases (dashed if between grid lines), shaded = scale strips (darker where they cross). |
| `lines.csv` | Every line on each axis: grid index, fraction of the paper, mm from each edge, and its kind (`edge`, `boundary`, `scale`, or `scale-sub` for a scale crease between grid lines). |
| `summary.txt` | Paper and scale sizes, each chunk's size and share compared with the target, and how to divide the square into N. |

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
