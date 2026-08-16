# Persistent Local Cohomology of a Point Cloud

Reference implementation of **He–Wei, *Persistent Local Cohomology of
Stanley–Reisner Rings at Vertex Primes***.

**In one sentence:** for a point cloud and a chosen point *i*, the program
computes a barcode that describes the shape of the cloud *around point i*
at every scale — detecting outliers, local dimension, branch points,
boundary points, and cut points, one point at a time.

## Requirements

- Python 3.8+ — **no dependencies** for the computation itself
  (pure standard library).
- `matplotlib` — only if you want the PNG plots (`pip install matplotlib`;
  otherwise pass `--no-plot`).

## Quick start (30 seconds)

```bash
sh reproduce.sh
```

runs the full test suite, reproduces every explicit number in the paper,
runs the two demos and the real-data example, and checks the outputs
against `expected_output/`. Or run one example by hand:

```bash
python3 run_plc.py --cloud demo/9-cut-points.csv --points connector --maxdim 8 --out my_run
```

## Input and output

The algorithm takes **two arguments**:

| argument | what | how to pass it |
|---|---|---|
| 1. the point cloud | coordinates (`--cloud points.csv`, one point per row, optional label column) **or** a pairwise-distance matrix (`--matrix D.csv`, optional labels) | see `demo/9-cut-points.csv` and `examples/benzene_D.csv` |
| 2. the points to analyze | `--points connector,A2` (labels or indices) | default: all points |

Output per run:

- `<out>_barcodes.csv` — one row per bar:
  `vertex, summand, sigma, q, birth, death, essential`
- `<out>_barcodes.png` — barcode plot (unless `--no-plot`)

## How to read a barcode (the 60-second version)

For each analyzed point *i*, the cloud is grown into a family of simplicial
complexes (Vietoris–Rips: connect points at distance ≤ *t*, for every *t*).
The invariant tracks two complexes attached to *i*:

- **link** — what the *neighbors* of *i* look like;
- **deletion (del)** — what the neighborhood looks like *with i removed*.

Each **bar** `[birth, death)` is one topological feature of the link/del
complex with its lifespan in distance units; the **color** `q` is the
degree (q=1: an extra connected piece; q=2: a loop; q=3: a cavity; feature
dimension is q−1). The number of bars of color `q` alive at scale `t` is
exactly a vector-space dimension: dim H̃_{q−1} of the link/del complex at
scale `t` — equivalently the dimension of the corresponding graded piece
of the local cohomology module `H^q_{p_i}` at the vertex prime
`p_i = (x_j : j ≠ i)`. A bar marked **essential** (`death = inf`) survives
into the final complex.

What the bars detect, per point:

| signature | meaning |
|---|---|
| long link bar, q = 0 | isolated point / outlier (length = isolation score) |
| dominant link bar at q = d | point interior to a d-dimensional region |
| m extra link bars at q = 1 | branch point where m+1 filaments meet |
| long **del** bar at q = 1 | **cut point**: removing it disconnects its neighborhood |
| link bars no manifold point can have | singular point (e.g. a cone apex — the paper's example) |

The `sigma` column refines the invariant: `{}` is the coarse slice
(the full link/del of *i*); non-empty faces probe the structure transverse
to a face at *i* (`--sigma all` computes them; default is `{}` only).

## The included examples

1. **`demo/`** — a 9-point cloud: two 4-point clusters bridged by one
   `connector` point. The demo detects the connector as a cut point: with
   it removed, its neighborhood stays disconnected over the scale window
   [0, 5) (deletion bar), and its link splits in two over [2.55, 5)
   (link bar). A comparison run on an ordinary cluster point shows no such
   bars. Full walkthrough: `demo/instruction.md`.
2. **`examples/benzene_D.csv`** — real data: the carbon–carbon distance
   matrix of benzene (Å). Run with `--two-prime` to also certify the
   result is characteristic-independent.
3. **`--demo`** — the paper's cone example (a non-flag filtration, given
   as an explicit filtered complex), reproducing the paper's Figure-4
   table bar for bar, including the nonempty-sigma essential bar.

## Options you may want

| flag | meaning |
|---|---|
| `--sigma empty\|all` | coarse slice only (default) / full face-indexed collection |
| `--maxdim D` | cap simplex dimension (default 3). Spurious *essential* bars can appear exactly at degrees q = D (link) and q = D+1 (del); the program warns about them. Raise `--maxdim` to check a suspect bar. |
| `--max-edge T` | stop the filtration at distance T |
| `--two-prime` | compute over two large primes; agreement certifies no torsion |
| `--exact` | independent exact (rational) backend, for small inputs |
| `--self-test` | reproduce every explicit bar in the paper and exit |


## Files

| file | role |
|---|---|
| `persistent_local_cohomology.py` | the library (pure Python, no dependencies) |
| `run_plc.py` | command-line interface |
| `test_plc.py` | acceptance tests |
| `reproduce.sh` | one-command full reproduction |
| `demo/` | cut-point demo: input + step-by-step instructions |
| `examples/` | benzene distance matrix (real data) |
| `expected_output/` | reference CSVs that `reproduce.sh` checks against |

## Citation

If you use this code, please cite the paper (He–Suwayyid-Wei, *Persistent Local
Cohomology of Stanley–Reisner Rings at Vertex Primes*).

## License

MIT — see `LICENSE`.
