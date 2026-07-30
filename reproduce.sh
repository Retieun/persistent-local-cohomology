#!/bin/sh
# Reproduce every result in this folder from scratch.
# Requires: Python 3.8+; matplotlib (optional, for the PNG plots).
set -e
cd "$(dirname "$0")"

echo "== 1/5  acceptance tests (12 checks, incl. fast==exact and stability) =="
python3 test_plc.py

echo
echo "== 2/5  paper self-test: every explicit bar in the paper =="
python3 run_plc.py --self-test

echo
echo "== 3/5  paper cone example, full sigma collection =="
python3 run_plc.py --demo --sigma all --out out_cone --no-plot

echo
echo "== 4/5  cut-point demo (9-cut-points.csv, connector) =="
python3 run_plc.py --cloud demo/9-cut-points.csv --points connector --maxdim 8 --out out_cut_point --no-plot
python3 run_plc.py --cloud demo/9-cut-points.csv --points A2 --maxdim 8 --out out_cluster_A2 --no-plot

echo
echo "== 5/5  real data: benzene (C-C distance matrix, Angstroms) =="
python3 run_plc.py --matrix examples/benzene_D.csv --two-prime --out out_benzene --no-plot

echo
echo "== comparing against expected_output/ =="
python3 - <<'EOF'
import csv, math, sys

def load(p):
    rows = sorted(tuple(r) for r in csv.reader(open(p)))
    return rows

pairs = [("out_cut_point_barcodes.csv",  "expected_output/cut_point_demo_barcodes.csv"),
         ("out_cluster_A2_barcodes.csv", "expected_output/cluster_point_A2_barcodes.csv"),
         ("out_cone_barcodes.csv",       "expected_output/cone_barcodes.csv")]
ok = True
for got, exp in pairs:
    same = load(got) == load(exp)
    print(("MATCH   " if same else "MISMATCH") + f"  {got}  vs  {exp}")
    ok = ok and same
sys.exit(0 if ok else 1)
EOF

echo
echo "All reproduction steps finished. Outputs: out_*.csv"
