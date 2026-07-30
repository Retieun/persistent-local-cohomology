#!/usr/bin/env python3
"""
Persistent Local Cohomology of Stanley-Reisner Rings at Vertex Primes
=====================================================================
Reference implementation of He-Wei, *Persistent Local Cohomology of
Stanley-Reisner Rings at Vertex Primes* (draft3.tex).  Implements the full
grand plan (GRAND_PLAN.md): general-sigma barcodes, explicit filtered-complex
input, lambda-numbers, exact cross-check backend, and stability utilities.

THE INVARIANT
    For a monotone filtration Delta^t = {tau : f(tau) <= t} of a finite
    simplicial complex on [n], and a vertex i, the persistent local cohomology
    of k[Delta^.] at the vertex prime p_i = (x_j : j != i) decomposes
    (Thm thm:structure + Prop prop:reduction) into a finite face-indexed
    collection of one-parameter persistence modules, indexed by
        (summand, sigma)  with  summand in {del, link},
        sigma in del_Delta(i)   for the deletion summand  (x_i-degree c = 0),
        sigma in link_Delta(i)  for the link summand      (x_i-degree c >= 1).
    Each module is interval-decomposable (Thm thm:finite-type); its barcode
    B^q_{i,del,sigma} / B^q_{i,link,sigma} (Def def:barcode) is computed here.

THE ALGORITHM (persistent links-Hochster formula, Thm thm:persistent-hochster)
    With d = q - |sigma| - 1, the barcode equals the persistent REDUCED
    homology in degree d of a derived filtered complex:
        deletion:  K_del^t  = lk_{Delta^t}(sigma) ∩ del_{Delta^t}(i),
                   a face tau enters at f(tau ∪ sigma);
        link:      K_link^t = lk_{Delta^t}(sigma ∪ {i}),
                   a face tau enters at f(tau ∪ sigma ∪ {i}).
    The enumeration is the bijection F |-> F \\ sigma (resp. F \\ (sigma∪{i}))
    over faces F of Delta containing sigma (resp. sigma ∪ {i}); the empty
    face of the derived complex therefore enters at f(sigma) (resp.
    f(sigma ∪ {i})) -- the void-complex convention of the theorem: nothing
    exists before the indexing face itself appears.

COEFFICIENTS
    Fast backend: single-pass persistence reduction over F_p (large prime,
    default 2^31 - 1).  This IS the paper's invariant for k = F_p.  Boundary
    matrices are 0/±1, so two large primes disagree exactly when p-torsion is
    present; agreement certifies the characteristic-0 barcode (two_prime).
    Exact backend: per-pair induced ranks over Q (pure-Fraction linear
    algebra, no dependencies) -- an algorithmically independent cross-check.

Distance units of the input carry through to all barcode endpoints.
"""
from __future__ import annotations
from itertools import combinations
from fractions import Fraction
import csv
import math

DEFAULT_PRIME = 2147483647           # 2^31 - 1
SECOND_PRIME = 2305843009213693951   # 2^61 - 1
INF = float('inf')


# ==========================================================================
# FilteredComplex
# ==========================================================================
class FilteredComplex:
    """A finite simplicial complex Delta with a monotone f: Delta -> R.

    Stores faces as frozensets of integer vertex ids, f-values as floats.
    The empty face is always present; its f-value defaults to the minimum
    vertex value (the filtration starts when the first vertex appears).
    """

    def __init__(self, pairs, labels=None, validate=True):
        faces = {}
        for face, val in pairs:
            fs = frozenset(int(v) for v in face)
            v = float(val)
            if fs in faces and faces[fs] != v:
                raise ValueError(f"face {sorted(fs)} given twice with "
                                 f"different values {faces[fs]} and {v}")
            faces[fs] = v
        verts = sorted({v for f in faces for v in f})
        if frozenset() not in faces:
            faces[frozenset()] = min(
                (faces[frozenset([u])] for u in verts
                 if frozenset([u]) in faces), default=0.0)
        self.faces = faces
        self.vertices = verts
        if labels is None:
            labels = [str(v) for v in verts]
        self.labels = {v: str(l) for v, l in zip(verts, labels)}
        if validate:
            self._validate()

    def _validate(self):
        for f, val in self.faces.items():
            if not f:
                continue
            fs = sorted(f)
            for k in range(len(fs)):
                sub = frozenset(fs[:k] + fs[k + 1:])
                if sub not in self.faces:
                    raise ValueError(
                        f"not a simplicial complex: face {fs} present but "
                        f"subface {sorted(sub)} missing")
                if self.faces[sub] > val + 1e-12:
                    raise ValueError(
                        f"f is not monotone: f({sorted(sub)})="
                        f"{self.faces[sub]} > f({fs})={val}")

    @property
    def dim(self):
        return max((len(f) - 1 for f in self.faces), default=-1)

    def del_faces(self, i):
        """Faces of del_Delta(i) (all faces not containing i, incl. empty)."""
        return [f for f in self.faces if i not in f]

    def link_faces(self, i):
        """Faces of link_Delta(i) = {tau : i not in tau, tau ∪ {i} in Delta}."""
        return [f - {i} for f in self.faces if i in f]

    def index_of(self, token):
        """Resolve a vertex label or integer-string to a vertex id."""
        tok = str(token).strip()
        for v, l in self.labels.items():
            if l == tok:
                return v
        return int(tok)


def vietoris_rips(dmat=None, points=None, max_dim=None, max_edge=None,
                  labels=None):
    """Vietoris-Rips (flag) filtration: f(tau) = diam(tau); a simplex enters
    once all its edges are present.  `dmat` is a symmetric zero-diagonal
    distance matrix; alternatively give `points` (coordinates, Euclidean).
    `max_dim` caps simplex dimension (default min(n-1, 3)); `max_edge` stops
    the filtration at that distance.
    """
    if dmat is None:
        if points is None:
            raise ValueError("give dmat or points")
        pts = [list(map(float, p)) for p in points]
        n = len(pts)
        dmat = [[math.dist(pts[a], pts[b]) for b in range(n)]
                for a in range(n)]
    D = [list(map(float, row)) for row in dmat]
    n = len(D)
    validate_matrix(D)
    if max_dim is None:
        max_dim = min(n - 1, 3)
    pairs = []
    for size in range(1, max_dim + 2):
        for c in combinations(range(n), size):
            b = _diam(D, c)
            if max_edge is None or b <= max_edge:
                pairs.append((c, b))
    fc = FilteredComplex(pairs, labels=labels, validate=False)
    fc.rips_maxdim = max_dim          # remembered for truncation warnings
    return fc


def _diam(D, tup):
    if len(tup) < 2:
        return 0.0
    return max(D[a][b] for a, b in combinations(tup, 2))


def validate_matrix(D):
    """Raise ValueError on inputs that do not yield a usable filtration."""
    n = len(D)
    if any(len(row) != n for row in D):
        raise ValueError("matrix is not square")
    for a in range(n):
        if abs(D[a][a]) > 1e-12:
            raise ValueError(f"diagonal entry D[{a}][{a}]={D[a][a]} is "
                             "nonzero (expected zero-diagonal distances)")
        for b in range(a + 1, n):
            if abs(D[a][b] - D[b][a]) > 1e-9:
                raise ValueError(f"matrix not symmetric at ({a},{b}): "
                                 f"{D[a][b]} vs {D[b][a]}")
    offdiag = {D[a][b] for a in range(n) for b in range(a + 1, n)}
    if offdiag <= {0.0, 1.0}:
        raise ValueError(
            "matrix looks binary/unweighted (off-diagonal entries only 0/1). "
            "A single graph has no filtration parameter, so there is no "
            "persistence to compute. Provide real distances/weights, or lift "
            "the graph to its shortest-path metric first.")


# ==========================================================================
# Derived filtrations (general sigma)
# ==========================================================================
def K_del_filtration(fc, i, sigma=frozenset()):
    """Filtered faces of K_del = lk(sigma) ∩ del(i): pairs (tau, birth) with
    tau = F \\ sigma over faces F ⊇ sigma, i ∉ F; birth = f(F) = f(tau∪sigma).
    Includes the empty face at f(sigma) (void-complex convention)."""
    sigma = frozenset(sigma)
    if i in sigma:
        raise ValueError(f"sigma {sorted(sigma)} contains the vertex {i}")
    if sigma not in fc.faces:
        raise ValueError(f"sigma {sorted(sigma)} is not a face of Delta")
    return [(tuple(sorted(F - sigma)), val) for F, val in fc.faces.items()
            if i not in F and sigma <= F]


def K_link_filtration(fc, i, sigma=frozenset()):
    """Filtered faces of K_link = lk(sigma ∪ {i}): pairs (tau, birth) with
    tau = F \\ (sigma∪{i}) over faces F ⊇ sigma∪{i}; birth = f(F).
    Includes the empty face at f(sigma ∪ {i})."""
    si = frozenset(sigma) | {i}
    if si not in fc.faces:
        raise ValueError(f"sigma ∪ {{i}} = {sorted(si)} is not a face of Delta")
    return [(tuple(sorted(F - si)), val) for F, val in fc.faces.items()
            if si <= F]


# ==========================================================================
# Fast backend: single-pass persistence reduction over F_p
# ==========================================================================
def reduced_persistence(faces, p=DEFAULT_PRIME):
    """Persistent REDUCED homology of a filtered complex given as a list of
    (vertex-tuple, birth) pairs (must include the empty face () and be closed
    under taking subsets).  Standard column reduction over F_p.
    Returns [(homological_dim, birth, death_or_None)], zero-length bars
    dropped, unpaired columns reported as essential (death None)."""
    order = sorted(range(len(faces)),
                   key=lambda j: (faces[j][1], len(faces[j][0]), faces[j][0]))
    pos = {tuple(faces[j][0]): r for r, j in enumerate(order)}
    filt = [faces[j][1] for j in order]
    dim = [len(faces[j][0]) - 1 for j in order]
    cols = []
    for j in order:
        f = faces[j][0]
        col = {}
        if len(f) >= 1:
            for k in range(len(f)):
                sub = f[:k] + f[k + 1:]
                col[pos[sub]] = ((-1) ** k) % p
        cols.append(col)
    low_map = {}
    R = [None] * len(order)
    neg_low = {}
    for r in range(len(order)):
        col = dict(cols[r])
        while col:
            l = max(col)
            if l not in low_map:
                break
            k = low_map[l]
            factor = (col[l] * pow(R[k][l], p - 2, p)) % p    # Fermat inverse
            for row, cf in R[k].items():
                nv = (col.get(row, 0) - factor * cf) % p
                if nv:
                    col[row] = nv
                elif row in col:
                    del col[row]
        R[r] = col
        if col:
            l = max(col)
            low_map[l] = r
            neg_low[r] = l
    paired = set(neg_low.values())
    bars = []
    for r, l in neg_low.items():
        if filt[r] > filt[l]:
            bars.append((dim[l], filt[l], filt[r]))
    for r in range(len(order)):
        if not R[r] and r not in paired:
            bars.append((dim[r], filt[r], None))
    return bars


# ==========================================================================
# Exact backend: per-pair induced ranks over Q (Fraction linear algebra)
# ==========================================================================
def _rref(mat):
    """Row-reduce a list-of-lists of Fractions in place; return pivot cols."""
    if not mat:
        return []
    rows, ncols = len(mat), len(mat[0])
    pivots, r = [], 0
    for c in range(ncols):
        piv = next((k for k in range(r, rows) if mat[k][c] != 0), None)
        if piv is None:
            continue
        mat[r], mat[piv] = mat[piv], mat[r]
        inv = Fraction(1, 1) / mat[r][c]
        mat[r] = [x * inv for x in mat[r]]
        for k in range(rows):
            if k != r and mat[k][c] != 0:
                f = mat[k][c]
                mat[k] = [a - f * b for a, b in zip(mat[k], mat[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return pivots


def _rank_frac(cols):
    """Rank of a matrix given as a list of columns (lists of ints/Fractions)."""
    if not cols:
        return 0
    mat = [[Fraction(c[r]) for c in cols] for r in range(len(cols[0]))]
    return len(_rref(mat))


def _nullspace_frac(cols, nrows):
    """Nullspace basis of the matrix with the given columns (each of length
    nrows); vectors returned in the column-index coordinates."""
    ncols = len(cols)
    if ncols == 0:
        return []
    mat = [[Fraction(cols[c][r]) for c in range(ncols)] for r in range(nrows)]
    pivots = _rref(mat)
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc_ in free:
        v = [Fraction(0)] * ncols
        v[fc_] = Fraction(1)
        for r, pc in enumerate(pivots):
            v[pc] = -mat[r][fc_]
        basis.append(v)
    return basis


def _boundary_cols(faces_set, d, basis_d, basis_dm1):
    """Columns of the reduced boundary map C_d -> C_{d-1} for the complex
    whose faces include `basis_d` (d-faces) and `basis_dm1` ((d-1)-faces).
    d = 0 maps vertices to the empty face with coefficient 1."""
    ri = {f: k for k, f in enumerate(basis_dm1)}
    cols = []
    for f in basis_d:
        col = [0] * len(basis_dm1)
        fs = sorted(f)
        for k in range(len(fs)):
            sub = frozenset(fs[:k] + fs[k + 1:])
            col[ri[sub]] += (-1) ** k
        cols.append(col)
    return cols


def exact_reduced_persistence(faces):
    """Same contract as reduced_persistence, but over Q by per-pair induced
    ranks (algorithmically independent of the reduction pairing).  Intended
    for small cross-checks only: O(T^2) rank computations."""
    levels = sorted({b for _, b in faces})
    T = len(levels) - 1
    fsets = [frozenset(f) for f, _ in faces]
    births = [b for _, b in faces]
    maxd = max((len(f) - 1 for f, _ in faces), default=-1)
    # faces present at each level, grouped by dimension
    by_level = []
    for t in levels:
        present = {}
        for fs, b in zip(fsets, births):
            if b <= t:
                present.setdefault(len(fs) - 1, []).append(fs)
        for dd in present:
            present[dd].sort(key=lambda f: tuple(sorted(f)))
        by_level.append(present)
    bars = []
    for d in range(-1, maxd + 1):
        # cycle spaces Z_d(K^s) embedded in K^T's d-chain coordinates,
        # boundary images B_d(K^t) in the same coordinates
        top_basis = by_level[T].get(d, [])
        idx_top = {f: k for k, f in enumerate(top_basis)}
        Z, Bcols = [], []
        for s in range(T + 1):
            bd = by_level[s].get(d, [])
            bdm1 = by_level[s].get(d - 1, [])
            cols = _boundary_cols(None, d, bd, bdm1)
            ns = _nullspace_frac(cols, len(bdm1))
            emb = []
            for v in ns:
                w = [Fraction(0)] * len(top_basis)
                for c, f in enumerate(bd):
                    w[idx_top[f]] = v[c]
                emb.append(w)
            Z.append(emb)
            bd1 = by_level[s].get(d + 1, [])
            colsB = _boundary_cols(None, d + 1, bd1, bd)
            embB = []
            for col in colsB:
                w = [Fraction(0)] * len(top_basis)
                for r, f in enumerate(bd):
                    w[idx_top[f]] = col[r]
                embB.append(w)
            Bcols.append(embB)
        rankB = [_rank_frac(B) for B in Bcols]

        def b(s, t):
            if s < 0 or t < 0 or s > t:
                return 0
            return _rank_frac(Z[s] + Bcols[t]) - rankB[t]

        Bcache = {(s, t): b(s, t) for s in range(T + 1)
                  for t in range(s, T + 1)}

        def bb(s, t):
            return 0 if (s < 0 or t < 0 or s > t) else Bcache[(s, t)]

        for bi in range(T + 1):
            for di in range(bi + 1, T + 1):
                mu = ((bb(bi, di - 1) - bb(bi, di))
                      - (bb(bi - 1, di - 1) - bb(bi - 1, di)))
                for _ in range(max(mu, 0)):
                    bars.append((d, levels[bi], levels[di]))
            mu = bb(bi, T) - bb(bi - 1, T)
            for _ in range(max(mu, 0)):
                bars.append((d, levels[bi], None))
    return bars


# ==========================================================================
# Top level: per-vertex barcodes over the face-indexed collection
# ==========================================================================
def _sigma_range(fc, i, summand, faces_sigma, sigma_maxsize, max_q):
    """Admissible sigma for (i, summand), per prop:reduction:
    deletion: sigma in del_Delta(i);  link: sigma in link_Delta(i).
    Enumerated from the final complex Delta.  Since the minimal possible
    cohomological degree of a (i, summand, sigma) barcode is |sigma|
    (d >= -1  =>  q = d + |sigma| + 1 >= |sigma|), sigma with |sigma| > max_q
    are pruned."""
    if faces_sigma == "empty":
        cands = [frozenset()]
    elif faces_sigma == "all":
        cands = (fc.del_faces(i) if summand == 'del' else fc.link_faces(i))
    else:
        cands = [frozenset(s) for s in faces_sigma]
    out = []
    for s in cands:
        s = frozenset(s)
        if i in s:
            continue
        if summand == 'del' and s not in fc.faces:
            continue
        if summand == 'link' and (s | {i}) not in fc.faces:
            continue
        if sigma_maxsize is not None and len(s) > sigma_maxsize:
            continue
        if max_q is not None and len(s) > max_q:
            continue
        out.append(s)
    return sorted(set(out), key=lambda s: (len(s), tuple(sorted(s))))


def per_vertex_barcodes(fc, points=None, faces_sigma="empty", max_q=None,
                        backend="fast", p=DEFAULT_PRIME, sigma_maxsize=None):
    """Per-vertex local-cohomology barcodes of the filtration `fc`.

    Arguments
        fc          : FilteredComplex (from vietoris_rips(...) or explicit)
        points      : vertices to analyze (ids); default all
        faces_sigma : "empty" (sigma = ∅ only), "all" (the full face-indexed
                      collection of def:barcode), or an explicit iterable of
                      faces
        max_q       : cap on cohomological degree q (default dim(fc) + 1)
        backend     : "fast" (F_p reduction) or "exact" (Q, per-pair ranks)
        p           : prime for the fast backend
        sigma_maxsize : cap |sigma| in "all" mode

    Returns  res[(i, summand, sigma)] = {q: [(birth, death_or_None), ...]}
    with summand in {"del", "link"}, sigma a frozenset, bars sorted, and the
    degree shift q = d + |sigma| + 1 applied (thm:persistent-hochster).
    """
    if points is None:
        points = list(fc.vertices)
    if max_q is None:
        max_q = fc.dim + 1
    engine = (exact_reduced_persistence if backend == "exact"
              else lambda F: reduced_persistence(F, p))
    res = {}
    for i in points:
        for summand, K in (('del', K_del_filtration),
                           ('link', K_link_filtration)):
            for sigma in _sigma_range(fc, i, summand, faces_sigma,
                                      sigma_maxsize, max_q):
                shift = len(sigma) + 1
                rec = {}
                for (d, b, dth) in engine(K(fc, i, sigma)):
                    q = d + shift
                    if 0 <= q <= max_q:
                        rec.setdefault(q, []).append((b, dth))
                for q in rec:
                    rec[q].sort(key=lambda z: (z[0], INF if z[1] is None
                                               else z[1]))
                res[(i, summand, sigma)] = rec
    return res


def torsion_report(fc, points=None, faces_sigma="empty", max_q=None,
                   p1=DEFAULT_PRIME, p2=SECOND_PRIME, sigma_maxsize=None):
    """Compute over two primes; return (bc_p1, disagreements), each
    disagreement a key (i, summand, sigma, q) whose bars differ between the
    fields -> evidence of p-torsion.  Agreement certifies the
    characteristic-0 barcode (boundary matrices are 0/±1 integer matrices)."""
    a = per_vertex_barcodes(fc, points, faces_sigma, max_q, "fast", p1,
                            sigma_maxsize)
    b = per_vertex_barcodes(fc, points, faces_sigma, max_q, "fast", p2,
                            sigma_maxsize)
    dis = []
    for key in set(a) | set(b):
        ra, rb = a.get(key, {}), b.get(key, {})
        for q in set(ra) | set(rb):
            if ra.get(q, []) != rb.get(q, []):
                dis.append(key + (q,))
    return a, sorted(dis, key=str)


# ==========================================================================
# Lambda numbers (def:numbers) and static dimensions
# ==========================================================================
def persistence_number(bars, s, t):
    """lambda^{q,s->t} from the barcode of one (i, summand, sigma, q) slice:
    the number of bars [b, d) with b <= s and t < d  (def:barcode <->
    def:numbers dictionary; s <= t).  With s = t this is the static dimension
    dim_k H^q at level s (cor:static-comb)."""
    if s > t:
        raise ValueError("need s <= t")
    return sum(1 for (b, d) in bars if b <= s and (d is None or d > t))


# ==========================================================================
# Stability utilities (thm:stability, def:coll-metric)
# ==========================================================================
def _bipartite_feasible(c1, c2, cost, eps):
    """Perfect-matching feasibility for bottleneck: every bar may match a bar
    of the other diagram (cost <= eps) or the diagonal (half-persistence
    <= eps).  Simple augmenting-path matching; sizes here are small."""
    n1, n2 = len(c1), len(c2)
    half1 = [(d - b) / 2 for (b, d) in c1]
    half2 = [(d - b) / 2 for (b, d) in c2]
    # bars of c2 that must be matched to a real bar of c1
    need2 = [j for j in range(n2) if half2[j] > eps]
    adj = {j: [k for k in range(n1) if cost[k][j] <= eps] for j in need2}
    match1 = [None] * n1

    def augment(j, seen):
        for k in adj[j]:
            if k in seen:
                continue
            seen.add(k)
            if match1[k] is None or augment(match1[k], seen):
                match1[k] = j
                return True
        return False

    matched2 = set()
    for j in need2:
        if not augment(j, set()):
            return False
        matched2.add(j)
    # bars of c1 unmatched so far must be diagonal-killable or matchable to
    # remaining c2 bars; run symmetric pass for the ones that need it
    used2 = {match1[k] for k in range(n1) if match1[k] is not None}
    for k in range(n1):
        if match1[k] is None and half1[k] > eps:
            found = False
            for j in range(n2):
                if j not in used2 and cost[k][j] <= eps:
                    used2.add(j)
                    found = True
                    break
            if not found:
                return False
    return True


def bottleneck_distance(bars1, bars2):
    """Bottleneck distance between two barcodes (lists of (birth, death),
    death None = essential).  Essential bars must match essential bars
    (unequal counts -> inf); on the line their optimal minimax matching is
    the sorted one.  Finite bars: exact combinatorial bottleneck (binary
    search over candidate costs + matching feasibility)."""
    ess1 = sorted(b for (b, d) in bars1 if d is None)
    ess2 = sorted(b for (b, d) in bars2 if d is None)
    if len(ess1) != len(ess2):
        return INF
    ess_cost = max((abs(a - b) for a, b in zip(ess1, ess2)), default=0.0)
    fin1 = [(b, d) for (b, d) in bars1 if d is not None]
    fin2 = [(b, d) for (b, d) in bars2 if d is not None]
    if not fin1 and not fin2:
        return ess_cost
    cost = [[max(abs(b1 - b2), abs(d1 - d2)) for (b2, d2) in fin2]
            for (b1, d1) in fin1]
    cands = {0.0, ess_cost}
    cands.update((d - b) / 2 for (b, d) in fin1 + fin2)
    cands.update(cost[k][j] for k in range(len(fin1))
                 for j in range(len(fin2)))
    lo, hi = 0.0, max(cands)
    feasible = sorted(c for c in cands if c >= 0)
    best = None
    for c in feasible:                       # sizes are small: linear scan
        if _bipartite_feasible(fin1, fin2, cost, c):
            best = c
            break
    return max(ess_cost, best if best is not None else INF)


def collection_distance(res_f, res_g, i, q):
    """d_coll of def:coll-metric for vertex i and degree q: the max over the
    face index set (summand, sigma) of the bottleneck distances between the
    q-slices.  res_f, res_g are per_vertex_barcodes outputs computed on the
    SAME complex Delta (the index set is fixed by Delta and i)."""
    keys = {k for k in set(res_f) | set(res_g) if k[0] == i}
    dmax, per = 0.0, {}
    for k in sorted(keys, key=str):
        bf = res_f.get(k, {}).get(q, [])
        bg = res_g.get(k, {}).get(q, [])
        d = bottleneck_distance(bf, bg)
        per[k] = d
        dmax = max(dmax, d)
    return dmax, per


# ==========================================================================
# Static Betti numbers (independent check for cor:static-comb)
# ==========================================================================
def static_reduced_betti(faces, t):
    """Reduced Betti numbers over Q of the derived complex at level t,
    computed directly from boundary-matrix ranks (independent of the
    persistence pairing).  Returns {d: betti_d}; void complex -> {}."""
    present = {}
    for f, b in faces:
        if b <= t:
            present.setdefault(len(f) - 1, []).append(frozenset(f))
    if not present:
        return {}
    for d in present:
        present[d].sort(key=lambda f: tuple(sorted(f)))
    out = {}
    maxd = max(present)
    for d in range(-1, maxd + 1):
        bd = present.get(d, [])
        if not bd:
            continue
        bdm1 = present.get(d - 1, [])
        bd1 = present.get(d + 1, [])
        rk_d = _rank_frac(_boundary_cols(None, d, bd, bdm1)) if bdm1 else 0
        rk_d1 = _rank_frac(_boundary_cols(None, d + 1, bd1, bd)) if bd1 else 0
        beta = len(bd) - rk_d - rk_d1
        if beta:
            out[d] = beta
    return out


# ==========================================================================
# Truncation-artifact warnings (Rips mode)
# ==========================================================================
def cap_artifact_warnings(res, fc):
    """A capped m-skeleton always has a spurious essential top-degree class
    (its killers live one dimension above the cap).  For every sigma these
    artifacts sit at the SAME cohomological boundary degrees
    (deletion q = maxdim + 1, link q = maxdim, since the degree shift
    q = d + |sigma| + 1 exactly compensates the reduced top dimension).
    Only meaningful for Rips-truncated complexes; silent when
    maxdim >= n - 1 (untruncated) or for explicit filtered-complex input."""
    warns = []
    maxdim = getattr(fc, 'rips_maxdim', None)
    if maxdim is None or maxdim >= len(fc.vertices) - 1:
        return warns
    for (i, summand, sigma), rec in res.items():
        boundary_q = maxdim + 1 if summand == 'del' else maxdim
        if any(d is None for (b, d) in rec.get(boundary_q, [])):
            warns.append(
                f"vertex {fc.labels[i]}, sigma={format_sigma(sigma, fc)}: "
                f"essential {summand} bar at q={boundary_q} is at the "
                f"maxdim={maxdim} truncation boundary -> likely a cap "
                f"artifact, not a real feature (raise --maxdim to check)")
    return warns


# ==========================================================================
# I/O helpers
# ==========================================================================
def format_sigma(sigma, fc=None):
    if not sigma:
        return "{}"
    names = ([fc.labels[v] for v in sorted(sigma)] if fc is not None
             else [str(v) for v in sorted(sigma)])
    return "{" + ",".join(names) + "}"


def load_matrix(path):
    """CSV distance matrix with optional label header row / label column."""
    rows = list(csv.reader(open(path)))

    def isnum(x):
        try:
            float(x)
            return True
        except ValueError:
            return False

    header = rows[0] and not all(isnum(x) for x in rows[0])
    if header:
        labels = rows[0][1:] if not isnum(rows[0][0]) else rows[0]
        body = rows[1:]
    else:
        labels = None
        body = rows
    rowlab = body and not isnum(body[0][0])
    M = []
    for r in body:
        cells = r[1:] if rowlab else r
        M.append([float(x) for x in cells])
    if labels is None:
        labels = [str(k) for k in range(len(M))]
    return M, [l.strip() for l in labels]


def load_cloud(path):
    """CSV point cloud: one point per row, coordinates in columns.

    Optional header row and optional label column (first column
    non-numeric), auto-detected -- the coordinate analogue of
    ``load_matrix``.  Returns (points, labels).
    """
    rows = [r for r in csv.reader(open(path)) if r]

    def isnum(x):
        try:
            float(x)
            return True
        except ValueError:
            return False

    if rows and not all(isnum(x) for x in rows[0]):
        body = rows[1:] if any(isnum(x) for x in rows[1]) else rows[1:]
    else:
        body = rows
    if not body:
        raise ValueError(f'{path}: no data rows')
    rowlab = not isnum(body[0][0])
    pts, labels = [], []
    for k, r in enumerate(body):
        cells = r[1:] if rowlab else r
        pts.append([float(x) for x in cells])
        labels.append(r[0].strip() if rowlab else str(k))
    dims = {len(p) for p in pts}
    if len(dims) != 1:
        raise ValueError(f'{path}: rows have inconsistent dimension {dims}')
    return pts, labels


def write_csv(path, res, fc):
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['vertex', 'summand', 'sigma', 'q', 'birth', 'death',
                    'essential'])
        for (i, summand, sigma) in sorted(res, key=lambda k: (
                k[0], k[1], len(k[2]), tuple(sorted(k[2])))):
            rec = res[(i, summand, sigma)]
            for q in sorted(rec):
                for (b, d) in rec[q]:
                    w.writerow([fc.labels[i], summand,
                                format_sigma(sigma, fc), q, f'{b:g}',
                                'inf' if d is None else f'{d:g}',
                                'yes' if d is None else 'no'])


def summarize(res, fc):
    lines = []
    byv = {}
    for (i, summand, sigma), rec in res.items():
        for q, bars in rec.items():
            if any(d is None for (b, d) in bars):
                byv.setdefault(i, []).append(
                    (summand, format_sigma(sigma, fc), q))
    for i in sorted({k[0] for k in res}):
        ess = sorted(byv.get(i, []))
        tag = ('; '.join(f'{s} sigma={sg} q={q}' for s, sg, q in ess)
               if ess else 'none')
        lines.append(f'  vertex {fc.labels[i]}: essential bars -> {tag}')
    return '\n'.join(lines)


# ==========================================================================
# Plotting (sigma = ∅ slice per analyzed vertex)
# ==========================================================================
def plot_barcodes(res, fc, path, title=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    qcol = {0: '#4C72B0', 1: '#DD8452', 2: '#C44E52', 3: '#8172B3',
            4: '#937860', 5: '#64B5CD'}
    pts = sorted({k[0] for k in res})
    empty = frozenset()
    tvals = [v for (b, d) in
             (bar for k in res if k[2] == empty
              for bars in res[k].values() for bar in bars)
             for v in ((b,) if d is None else (b, d))]
    tmax = max(tvals) if tvals else 1.0
    tmax = tmax if tmax > 0 else 1.0
    pad = 0.06 * tmax
    fig, axes = plt.subplots(1, len(pts), figsize=(3.0 * len(pts), 3.6),
                             sharex=True, squeeze=False)
    axes = axes[0]
    for ax, i in zip(axes, pts):
        y = 0
        yt, yl, band = [], [], []
        for gi, summand in enumerate(('del', 'link')):
            rec = res.get((i, summand, empty), {})
            g0 = y
            for q in sorted(rec):
                for (b, d) in rec[q]:
                    ess = d is None
                    xend = tmax + pad if ess else d
                    ax.plot([b, xend], [y, y], color=qcol.get(q, '0.3'),
                            lw=2.6 if ess else 2.0, solid_capstyle='round',
                            zorder=3)
                    if ess:
                        ax.annotate('', xy=(tmax + 2.2 * pad, y),
                                    xytext=(tmax + pad, y),
                                    arrowprops=dict(arrowstyle='-|>',
                                                    color=qcol.get(q, '0.3'),
                                                    lw=2.4))
                    y += 1
            if y > g0:
                yt.append((g0 + y - 1) / 2)
                yl.append('deletion' if summand == 'del' else 'link')
                band.append((g0 - 0.5, y - 0.5, gi))
            y += 0.7
        for (y0, y1, gi) in band:
            if gi % 2 == 0:
                ax.axhspan(y0, y1, color='0.94', zorder=0)
        ax.set_yticks(yt)
        ax.set_yticklabels(yl, fontsize=8)
        ax.set_title(f'vertex {fc.labels[i]}', fontsize=9)
        ax.set_xlabel('filtration distance $t$', fontsize=8)
        ax.set_xlim(-0.03 * tmax, tmax + 3 * pad)
        ax.set_ylim(-0.7, max(y, 1) - 0.3)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    qs = sorted({q for k in res if k[2] == empty for q in res[k]})
    leg = [Line2D([0], [0], color=qcol.get(q, '0.3'), lw=2.4,
                  label=f'$q={q}$') for q in qs]
    leg.append(Line2D([0], [0], color='0.3', lw=2.6, marker='>',
                      markersize=5, label=r'essential ($\to\infty$)'))
    axes[-1].legend(handles=leg, loc='upper right', frameon=False,
                    fontsize=7)
    fig.suptitle(title or 'Persistent per-vertex local-cohomology barcodes '
                 r'(face-$\varnothing$ slice)', fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return path


# ==========================================================================
# The paper's cone example (ex:cone / Fig. 4) -- a NON-FLAG filtration
# ==========================================================================
def cone_example():
    """The filtration of ex:cone: cone with apex 0 over the boundary of the
    triangle {1,2,3}.  Rim edges at t=1, apex edges at t=2, side triangles at
    t=3; the rim triangle {1,2,3} is never a face (non-flag!), keeping Delta
    a disk."""
    faces = {(0,): 0, (1,): 0, (2,): 0, (3,): 0,
             (1, 2): 1, (2, 3): 1, (1, 3): 1,
             (0, 1): 2, (0, 2): 2, (0, 3): 2,
             (0, 1, 2): 3, (0, 2, 3): 3, (0, 1, 3): 3}
    return FilteredComplex(faces.items())


#: Every explicit bar of the paper's Fig. 4 (sigma = ∅ slice, vertices 0, 1).
PAPER_CONE_BARS = {
    (0, 'link', 0): [(0.0, 2.0)],
    (0, 'link', 1): [(2.0, 3.0), (2.0, 3.0)],
    (0, 'link', 2): [(3.0, None)],
    (0, 'del', 1): [(0.0, 1.0), (0.0, 1.0)],
    (0, 'del', 2): [(1.0, None)],
    (1, 'link', 0): [(0.0, 1.0)],
    (1, 'link', 1): [(1.0, 3.0), (2.0, 3.0)],
    (1, 'del', 1): [(0.0, 1.0), (0.0, 2.0)],
    (1, 'del', 2): [(2.0, 3.0)],
}


def reproduce_paper_example(backend="fast", p=DEFAULT_PRIME):
    """Check the computed cone barcodes against every explicit bar in the
    paper: the sigma = ∅ table of ex:cone (Fig. 4) AND the nonempty-face
    assertion B^2_{1,link,{0,2}} = {[3, inf)}.  Returns (ok, details)."""
    fc = cone_example()
    res = per_vertex_barcodes(fc, points=[0, 1], faces_sigma="all",
                              backend=backend, p=p)
    details, ok = [], True
    empty = frozenset()
    for i in (0, 1):
        for summand in ('del', 'link'):
            got = res.get((i, summand, empty), {})
            qs = set(got) | {q for (ii, ss, q) in PAPER_CONE_BARS
                             if ii == i and ss == summand}
            for q in sorted(qs):
                g = got.get(q, [])
                e = sorted(PAPER_CONE_BARS.get((i, summand, q), []),
                           key=lambda z: (z[0], INF if z[1] is None else z[1]))
                if g != e:
                    ok = False
                    details.append(f'MISMATCH ({i},{summand},sigma={{}},q={q}): '
                                   f'got {g}, paper {e}')
    key = (1, 'link', frozenset({0, 2}))
    got2 = res.get(key, {}).get(2, [])
    if got2 != [(3.0, None)]:
        ok = False
        details.append(f'MISMATCH B^2_(1,link,{{0,2}}): got {got2}, '
                       f'paper [(3.0, None)]')
    return ok, details
