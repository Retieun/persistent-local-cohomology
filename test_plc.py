#!/usr/bin/env python3
"""
Acceptance tests for persistent_local_cohomology (GRAND_PLAN.md section 6).
Run:  python3 test_plc.py        (plain; no pytest dependency needed)
"""
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import persistent_local_cohomology as plc

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


def bars_eq(a, b):
    key = lambda z: (z[0], math.inf if z[1] is None else z[1])
    return sorted(a, key=key) == sorted(b, key=key)


def random_rips(n, seed, max_dim=None):
    rng = random.Random(seed)
    pts = [(rng.uniform(0, 2), rng.uniform(0, 2)) for _ in range(n)]
    return plc.vietoris_rips(points=pts, max_dim=max_dim)


# 1. Paper cone example, sigma = ∅ table (Fig. 4), fast backend ------------
@check("1. Fig.4 cone table, sigma = empty (fast backend)")
def _():
    ok, details = plc.reproduce_paper_example()
    assert ok, "\n".join(details)


# 2. Paper's nonempty-sigma assertion + exact backend on the cone ----------
@check("2. cone, sigma != empty: B^2_(1,link,{0,2}) = {[3,inf)}; "
       "exact backend agrees on the full collection")
def _():
    fc = plc.cone_example()
    fast = plc.per_vertex_barcodes(fc, [0, 1], "all", backend="fast")
    key = (1, 'link', frozenset({0, 2}))
    assert fast[key].get(2) == [(3.0, None)], fast.get(key)
    exact = plc.per_vertex_barcodes(fc, [0, 1], "all", backend="exact")
    assert set(fast) == set(exact)
    for k in fast:
        qs = set(fast[k]) | set(exact[k])
        for q in qs:
            assert bars_eq(fast[k].get(q, []), exact[k].get(q, [])), \
                (k, q, fast[k].get(q), exact[k].get(q))


# 3. fast == exact on random Rips inputs, sigma = all ----------------------
@check("3. fast == exact, random Rips, full sigma collection")
def _():
    for seed in (1, 2, 3):
        fc = random_rips(5, seed)
        fast = plc.per_vertex_barcodes(fc, faces_sigma="all", backend="fast")
        exact = plc.per_vertex_barcodes(fc, faces_sigma="all",
                                        backend="exact")
        assert set(fast) == set(exact)
        for k in fast:
            for q in set(fast[k]) | set(exact[k]):
                assert bars_eq(fast[k].get(q, []), exact[k].get(q, [])), \
                    (seed, k, q)


# 4. two-prime certificate exercises cleanly -------------------------------
@check("4. two-prime agreement on random data and on the cone")
def _():
    for fc in (random_rips(6, 7), plc.cone_example()):
        _, dis = plc.torsion_report(fc, faces_sigma="all")
        assert dis == [], dis


# 5. void-complex convention: no bar born before f(sigma) ------------------
@check("5. void-complex convention (birth >= f(sigma) / f(sigma+i))")
def _():
    for seed in (11, 12):
        fc = random_rips(6, seed)
        res = plc.per_vertex_barcodes(fc, faces_sigma="all", backend="fast")
        for (i, summand, sigma), rec in res.items():
            T = (fc.faces[sigma] if summand == 'del'
                 else fc.faces[sigma | {i}])
            for q, bars in rec.items():
                for (b, d) in bars:
                    assert b >= T - 1e-12, (i, summand, sigma, q, b, T)


# 6. static sanity: lambda(s,s) == direct Betti of the derived complex -----
@check("6. static dimensions match direct Betti numbers (cor:static-comb)")
def _():
    for seed in (21, 22):
        fc = random_rips(5, seed)
        res = plc.per_vertex_barcodes(fc, faces_sigma="all", backend="fast")
        levels = sorted({v for v in fc.faces.values()})
        for (i, summand, sigma), rec in res.items():
            K = (plc.K_del_filtration(fc, i, sigma) if summand == 'del'
                 else plc.K_link_filtration(fc, i, sigma))
            shift = len(sigma) + 1
            for t in levels[::2]:
                betti = plc.static_reduced_betti(K, t)
                for d, beta in betti.items():
                    q = d + shift
                    lam = plc.persistence_number(rec.get(q, []), t, t)
                    assert lam == beta, (i, summand, sigma, q, t, lam, beta)
                # and no phantom dimensions where betti is 0
                for q, bars in rec.items():
                    d = q - shift
                    if d not in betti:
                        lam = plc.persistence_number(bars, t, t)
                        assert lam == 0, (i, summand, sigma, q, t, lam)


# 7. truncation guard fires at capped maxdim, silent untruncated -----------
@check("7. cap-artifact warnings: fire when capped, silent at maxdim >= n-1")
def _():
    # 4 points, all pairwise close -> full tetrahedron if uncapped
    D = [[0, 1, 1, 1], [1, 0, 1, 1.1], [1, 1, 0, 1.2], [1, 1.1, 1.2, 0]]
    capped = plc.vietoris_rips(dmat=D, max_dim=1)
    res_c = plc.per_vertex_barcodes(capped)
    assert plc.cap_artifact_warnings(res_c, capped), \
        "expected warnings at maxdim=1"
    full = plc.vietoris_rips(dmat=D, max_dim=3)
    res_f = plc.per_vertex_barcodes(full)
    assert plc.cap_artifact_warnings(res_f, full) == []


# 8. lambda-numbers: bar-counting dictionary of def:numbers ----------------
@check("8. persistence_number: lambda^{q,s->t} = #{bars [b,d): b<=s, t<d}")
def _():
    bars = [(0.0, 1.0), (0.0, 2.0), (1.0, None)]
    assert plc.persistence_number(bars, 0.0, 0.0) == 2
    assert plc.persistence_number(bars, 0.0, 0.5) == 2
    assert plc.persistence_number(bars, 0.0, 1.0) == 1   # [0,1) dead at 1
    assert plc.persistence_number(bars, 1.0, 1.5) == 2
    assert plc.persistence_number(bars, 1.0, 2.5) == 1   # only the essential
    assert plc.persistence_number(bars, 0.5, 2.5) == 0   # essential born at 1


# 9. stability: d_coll <= ||f - g||_inf on perturbed filtrations -----------
@check("9. stability (thm:stability): d_coll <= ||f-g||_inf, numerically")
def _():
    rng = random.Random(31)
    fc = plc.cone_example()
    eps = 0.25
    pairs_g = [(tuple(sorted(f)), v + rng.uniform(0, eps))
               for f, v in fc.faces.items() if f]
    # enforce monotonicity of g by rounding up along inclusions
    g = dict(pairs_g)
    for f in sorted(g, key=len):
        for k in range(len(f)):
            sub = f[:k] + f[k + 1:]
            if sub in g:
                g[f] = max(g[f], g[sub])
    delta = max(abs(g[tuple(sorted(f))] - v)
                for f, v in fc.faces.items() if f)
    fcg = plc.FilteredComplex(g.items())
    rf = plc.per_vertex_barcodes(fc, [0, 1], "all")
    rg = plc.per_vertex_barcodes(fcg, [0, 1], "all")
    for i in (0, 1):
        for q in (0, 1, 2):
            dcoll, _per = plc.collection_distance(rf, rg, i, q)
            assert dcoll <= delta + 1e-9, (i, q, dcoll, delta)


# 10. bottleneck distance basic identities ---------------------------------
@check("10. bottleneck: identity 0, essential mismatch inf, known value")
def _():
    b1 = [(0.0, 2.0), (1.0, None)]
    assert plc.bottleneck_distance(b1, b1) == 0.0
    assert plc.bottleneck_distance(b1, [(0.0, 2.0)]) == math.inf
    # single bars, shifted: max(|db|,|dd|)
    assert abs(plc.bottleneck_distance([(0.0, 2.0)], [(0.5, 2.0)]) - 0.5) \
        < 1e-12
    # a tiny bar matches the diagonal at half-persistence
    assert abs(plc.bottleneck_distance([(1.0, 1.2)], []) - 0.1) < 1e-12


# 11. input validation ------------------------------------------------------
@check("11. validate_matrix rejects binary / asymmetric / nonzero-diagonal")
def _():
    for bad in ([[0, 1], [1, 0]],                       # binary
                [[0, 1.0], [2.0, 0]],                   # asymmetric
                [[0.5, 1.0], [1.0, 0]]):                # nonzero diagonal
        try:
            plc.validate_matrix([[float(x) for x in r] for r in bad])
            raise AssertionError(f"accepted {bad}")
        except ValueError:
            pass
    # non-monotone explicit filtration rejected
    try:
        plc.FilteredComplex([((0,), 1.0), ((1,), 1.0), ((0, 1), 0.5)])
        raise AssertionError("accepted non-monotone f")
    except ValueError:
        pass
    # non-closed complex rejected
    try:
        plc.FilteredComplex([((0,), 0.0), ((0, 1), 1.0)])
        raise AssertionError("accepted non-closed face set")
    except ValueError:
        pass


# 12. degree-shift bookkeeping: q = d + |sigma| + 1 -------------------------
@check("12. degree shift: minimal q of a sigma-slice is >= |sigma|")
def _():
    fc = random_rips(6, 41)
    res = plc.per_vertex_barcodes(fc, faces_sigma="all", backend="fast")
    for (i, summand, sigma), rec in res.items():
        for q in rec:
            assert q >= len(sigma), (i, summand, sigma, q)


def main():
    failed = 0
    for name, fn in CHECKS:
        try:
            fn()
            print(f'PASS  {name}')
        except AssertionError as e:
            failed += 1
            print(f'FAIL  {name}\n      {e}')
    print(f'\n{len(CHECKS) - failed}/{len(CHECKS)} checks passed')
    raise SystemExit(1 if failed else 0)


if __name__ == '__main__':
    main()
