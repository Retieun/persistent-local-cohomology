#!/usr/bin/env python3
"""
CLI for persistent per-vertex local cohomology (He-Wei, draft3.tex).

Inputs (argument 1: the point cloud / filtered complex)
    --matrix D.csv    CSV distance matrix, optional label row/column
                      (the benzene_D.csv format)
    input.json        JSON in any of three shapes:
                        {"faces":  [[[0], 0.0], [[0,1], 1.0], ...]}
                            explicit monotone filtered complex (may be
                            non-flag; the only way to express e.g. the
                            paper's cone example)
                        {"points": [[x,y,...], ...], "max_dim": 2,
                         "max_edge": 1.6}   -> Vietoris-Rips
                        {"dmat":   [[...], ...], "max_dim": 2,
                         "max_edge": 1.6}   -> Vietoris-Rips
                      optional "labels": [...]
    --demo            the paper's cone example (ex:cone / Fig. 4)

Argument 2: the points to analyze
    --points/--vertices  comma-separated labels or indices (default: all)

Outputs: <out>_barcodes.csv (all computed (vertex, summand, sigma, q) bars)
and <out>_barcodes.png (the sigma = ∅ slice figure).
"""
import argparse
import json
import sys

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import persistent_local_cohomology as plc


def build_complex(args):
    if args.demo:
        return plc.cone_example()
    if args.matrix:
        M, labels = plc.load_matrix(args.matrix)
        return plc.vietoris_rips(dmat=M, max_dim=args.maxdim,
                                 max_edge=args.max_edge, labels=labels)
    if args.cloud:
        pts, labels = plc.load_cloud(args.cloud)
        return plc.vietoris_rips(points=pts, max_dim=args.maxdim,
                                 max_edge=args.max_edge, labels=labels)
    if not args.input:
        raise SystemExit("give an input JSON, --matrix D.csv, or --demo "
                         "(-h for help)")
    spec = json.load(open(args.input))
    labels = spec.get("labels")
    if "faces" in spec:
        fc = plc.FilteredComplex(
            ((face, val) for face, val in spec["faces"]), labels=labels)
        return fc
    if "points" in spec:
        return plc.vietoris_rips(points=spec["points"],
                                 max_dim=args.maxdim or spec.get("max_dim"),
                                 max_edge=args.max_edge or spec.get("max_edge"),
                                 labels=labels)
    if "dmat" in spec:
        return plc.vietoris_rips(dmat=spec["dmat"],
                                 max_dim=args.maxdim or spec.get("max_dim"),
                                 max_edge=args.max_edge or spec.get("max_edge"),
                                 labels=labels)
    raise SystemExit('input JSON must contain "faces", "points", or "dmat"')


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Persistent per-vertex local cohomology of a '
                    'Stanley-Reisner ring (He-Wei).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input', nargs='?', help='input JSON (see module doc)')
    ap.add_argument('--matrix', help='CSV distance matrix')
    ap.add_argument('--cloud', help='CSV point cloud (one point per row, '
                    'optional label column; pairwise Euclidean distances '
                    'are computed)')
    ap.add_argument('--demo', action='store_true',
                    help="the paper's cone example (non-flag filtration)")
    ap.add_argument('--points', '--vertices', dest='points',
                    help='comma-separated labels or indices (default: all)')
    ap.add_argument('--sigma', default='empty', choices=['empty', 'all'],
                    help='face-index range: sigma = ∅ only, or the full '
                         'face-indexed collection (default: empty)')
    ap.add_argument('--sigma-maxsize', type=int, default=None,
                    help='cap |sigma| in --sigma all mode')
    ap.add_argument('--max-q', type=int, default=None,
                    help='cap cohomological degree q (default: dim+1)')
    ap.add_argument('--maxdim', type=int, default=None,
                    help='Rips: cap simplex dimension (default: min(n-1,3))')
    ap.add_argument('--max-edge', '--max-radius', dest='max_edge',
                    type=float, default=None,
                    help='Rips: stop the filtration at this distance')
    ap.add_argument('--exact', action='store_true',
                    help='exact rational backend (slow; small inputs only)')
    ap.add_argument('--p', '--prime', dest='p', type=int,
                    default=plc.DEFAULT_PRIME,
                    help='prime for the fast backend (default 2^31-1)')
    ap.add_argument('--two-prime', action='store_true',
                    help='also compute over 2^61-1 and report p-torsion')
    ap.add_argument('--out', default=None,
                    help='output prefix (default: plc, or cone for --demo)')
    ap.add_argument('--no-plot', action='store_true')
    ap.add_argument('--self-test', action='store_true',
                    help='reproduce every explicit bar of the paper '
                         '(Fig. 4 sigma=∅ table + the nonempty-sigma '
                         'assertion) and exit')
    a = ap.parse_args(argv)

    if a.self_test:
        ok, details = plc.reproduce_paper_example()
        print('paper example reproduced bar-for-bar:', ok)
        for d in details:
            print(' ', d)
        sys.exit(0 if ok else 1)

    fc = build_complex(a)
    if a.points:
        pts = [fc.index_of(t) for t in a.points.split(',')]
    else:
        pts = list(fc.vertices)
    out = a.out or ('cone' if a.demo else 'plc')

    backend = 'exact' if a.exact else 'fast'
    if a.two_prime and not a.exact:
        res, dis = plc.torsion_report(fc, pts, a.sigma, a.max_q,
                                      a.p, plc.SECOND_PRIME, a.sigma_maxsize)
        if dis:
            print('P-TORSION DETECTED (bars differ between primes) at:')
            for (i, s, sg, q) in dis:
                print(f'  vertex {fc.labels[i]}, {s}, '
                      f'sigma={plc.format_sigma(sg, fc)}, q={q}')
        else:
            print('two-prime check: no p-torsion; barcode equals the '
                  'characteristic-0 result')
    else:
        res = plc.per_vertex_barcodes(fc, pts, a.sigma, a.max_q, backend,
                                      a.p, a.sigma_maxsize)

    plc.write_csv(f'{out}_barcodes.csv', res, fc)
    wrote = [f'{out}_barcodes.csv']
    if not a.no_plot:
        plc.plot_barcodes(res, fc, f'{out}_barcodes.png')
        wrote.append(f'{out}_barcodes.png')
    print(f'backend         : {backend}')
    print(f'points analyzed : {[fc.labels[i] for i in pts]}')
    print(f'sigma range     : {a.sigma}')
    print('essential-bar summary:')
    print(plc.summarize(res, fc))
    for w in plc.cap_artifact_warnings(res, fc):
        print('  [warning]', w)
    print('wrote ' + ' and '.join(wrote))


if __name__ == '__main__':
    main()
