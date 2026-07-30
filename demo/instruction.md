# Demo: detecting a cut point with per-vertex persistent local cohomology

This folder is a self-contained example of the algorithm's input contract

```
(point_cloud.csv, point)   →   per-point local-cohomology barcodes
```

Here the pair is **(`9-cut-points.csv`, `connector`)**.

## The data

`9-cut-points.csv` — a point cloud with one point per row (`label, x, y`;
labels optional). Nine points in the plane: two unit-square clusters of
four points each, `A1..A4` centered at (−3, 0) and `B1..B4` at (+3, 0),
bridged by a single point `connector` at the origin:

```
A1  A2        B1  B2
        connector
A3  A4        B3  B4
```

Distances that matter: within a cluster ≈ 1–1.41; cluster → connector
≈ 2.55; the shortest edge **between** the clusters (not through the
connector) is 5. So for scales `t` in `[2.55, 5)` the two clusters are
joined *only through the connector* — it is a cut point of the Rips
complex in that range.

## Run it

From this folder (Python 3, no dependencies beyond matplotlib for the plot):

```bash
python3 ../run_plc.py --cloud 9-cut-points.csv --points connector --maxdim 8 --out cut_point_demo
```

- `--cloud 9-cut-points.csv` — argument 1, the point cloud (pairwise
  Euclidean distances are computed; a precomputed distance matrix via
  `--matrix` works identically).
- `--points connector` — argument 2, the point to analyze.
- `--maxdim 8` — with 9 points the full skeleton is cheap; this removes
  all truncation artifacts (omit it on larger clouds and heed the cap
  warnings instead).

Outputs: `cut_point_demo_barcodes.csv` (all bars) and
`cut_point_demo_barcodes.png` (barcode plot).

## What you should see

The two bars that certify the cut point (both in the σ = ∅ slice):

| bar | meaning |
|---|---|
| **deletion, q = 1: `[0, 5)`** | remove `connector` and its neighborhood falls into two pieces; the extra component lives until `t = 5`, when the clusters finally touch directly. `q = 1` deletion bars count components of `del(connector)` beyond the first (Čech-style, via `H̃⁰`, shifted to `q = d+1`). |
| **link, q = 1: `[2.55, 5)`** | from the moment `connector` acquires neighbors on both sides (`t ≈ 2.55`) until the clusters merge (`t = 5`), its link has **two** connected components — the A-side and the B-side. A point interior to a single blob never shows this. |

The bar *length* (here `5 − 0` and `5 − 2.55`) is the scale range over
which the point acts as a bridge — a quantitative "cut-point score".

## Compare with a non-cut point

```bash
python3 ../run_plc.py --cloud 9-cut-points.csv --points A2 --maxdim 8 --out cluster_point_A2
```

`A2` sits inside cluster A. Its longest deletion `q = 1` bar dies at
`2.55` (as soon as the connector re-joins what removing `A2` briefly
separates — nothing global), and its link `q = 1` bar is the tiny
`[1, 1.41)`. No bar survives on `[2.55, 5)`: removing `A2` never
disconnects anything at bridge scales. The contrast between the two CSV
files *is* the detection.

Reference outputs for both runs are included in `../expected_output/`
(`cut_point_demo_barcodes.csv`, `cluster_point_A2_barcodes.csv`).

## Notes

- Analyze several points at once: `--points connector,A2,B1`, or omit
  `--points` for all nine.
- Add `--two-prime` for the torsion certificate, `--exact` for the
  slower exact ℚ backend, `--sigma all` for the full face-indexed
  collection (see `../README.md`).
- The mathematics: the deletion bars are the `x_i`-degree-0 part and the
  link bars the `x_i`-degree ≥ 1 part of the persistent local cohomology
  `H^q_{p_i}(k[Δ^•])` at the vertex prime `p_i`; see `../README.md` and
  the paper (He–Wei, *Persistent Local Cohomology of Stanley–Reisner
  Rings at Vertex Primes*).
