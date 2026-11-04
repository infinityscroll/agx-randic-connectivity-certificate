#!/usr/bin/env python3
"""
verify_theorem.py -- exact comparisons of paths and proper double comets for R*a,
relevant to the product clause of Conjecture 4 in

  M. Aouchiche, P. Hansen, M. Zheng, "Variable neighborhood search for extremal graphs 19:
  Further conjectures and results about the Randic index",
  MATCH Commun. Math. Comput. Chem. 58 (2007) 83-102
  (restated as Conjecture 1.2 in X. Li, Y. Shi, L. Wang, MATCH Commun. Math. Comput. Chem.
  68 (2012) 843-849).

THEOREM.  Let n in {10, 11, 12, 13}.  For every double comet DC(k1,k2,p) on n vertices
(k1 >= 2, k2 >= 2, p >= 0, k1 + k2 + p + 2 = n)

        R(DC(k1,k2,p)) * a(DC(k1,k2,p))  >  R(P_n) * a(P_n).

For n = 14 the inequality is reversed for DC(2,2,8) (also certified below).

Notation.  R(G) = sum over edges uv of (d_u d_v)^(-1/2)  (Randic index);
a(G) = second-smallest eigenvalue of the Laplacian L = D - A (algebraic connectivity);
DC(k1,k2,p) = two centres c1, c2 joined by a path with p internal vertices, with k1 pendant
vertices attached to c1 and k2 pendant vertices attached to c2  (n = k1 + k2 + p + 2).

Final product comparisons use exact integer / rational arithmetic (int, fractions.Fraction).
Floating point supplies trial endpoints for Method A and a noncertifying closed-form
sanity check. A bad guess can make certification fail, never make it succeed wrongly.

  R : each term (d_u d_v)^(-1/2) is enclosed in a rational interval of width < 10^-40 using
      math.isqrt (exact when d_u d_v is a perfect square).

  a, Method A (Sylvester's law of inertia):  for rational t, the number of eigenvalues of L
      that are < t equals the number of negative pivots of the exact LDL^T factorisation of
      L - tI.  (No pivot can vanish when t is a non-integer rational: the leading principal
      minors of L - tI are +-monic integer polynomials in t, whose rational roots are integers.
      Zero pivots are nevertheless checked for and would abort.)  From a floating guess we
      pick decimal rationals lo < hi and certify  #{eig < lo} = 1  and  #{eig < hi} >= 2,
      i.e.  lo <= a < hi.

  a, Method B (Sturm sequences; no floating point at all):  the integer characteristic
      polynomial p(x) = det(xI - L) is computed exactly by the Faddeev-LeVerrier recursion
      over Z (checked via Cayley-Hamilton), q(x) = p(x)/x.  q(0) != 0 certifies that 0 is a
      simple eigenvalue, hence a = smallest root of q (all roots are real, L is symmetric).
      The smallest root is isolated by exact bisection of [0, n + 1/2] driven by Sturm
      sign-variation counts, then rounded outward to decimals and re-certified:
      no root of q in (-inf, lo]  and  at least one root of q in (lo, hi).

  P_n, Method C (closed form):  a(P_n) = 2 - 2 cos(pi/n), with rigorous rational bounds
      (Machin's formula for pi, alternating Taylor series for cos); used as a cross-check.

  Optional: if sympy is importable, its characteristic polynomial (Berkowitz) is compared
  with ours for every graph.  The certificate itself does not use sympy.

Usage:
    python3 verify_theorem.py              # theorem (n = 10..13) + certified reversal at n = 14
    python3 verify_theorem.py --n 14 15    # just the tables for the listed n
    python3 verify_theorem.py --no-sympy   # skip the optional sympy cross-check
Exit status 0 iff all requested comparisons are decided and all cross-checks pass.
The default invocation certifies the finite theorem and the reversal at n=14.
See docs/SCOPE.md for the definition issue in interpreting the original conjecture.
"""
import argparse
import platform
import sys
import time
from fractions import Fraction as Fr
from math import gcd, isqrt, sqrt

try:
    import sympy
except ImportError:              # sympy is optional
    sympy = None

R_DIGITS = 40     # each (d_u d_v)^(-1/2) is enclosed to within 10^-40
A_DIGITS = 14     # Method A test points lie on the grid 10^-14 Z (widened on failure)
B_DIGITS = 15     # Method B interval is rounded outward to the grid 10^-15 Z
SHOW = 13         # printed decimals: lower bounds rounded down, upper bounds rounded up


def require(condition, message):
    """A certificate check that remains active under python -O."""
    if not condition:
        raise RuntimeError(message)


# ======================================================================================
# Graphs
# ======================================================================================
def path_graph(n):
    """P_n on vertices 0..n-1."""
    return n, [(i, i + 1) for i in range(n - 1)]


def double_comet(k1, k2, p):
    """DC(k1,k2,p).  Vertex order: the k1 leaves of c1, c1, the p internal path vertices,
    the k2 leaves of c2, c2 (= last vertex).  Every vertex except the last then has exactly
    one later neighbour, so the LDL^T of Method A has no fill-in (this only affects speed;
    the eigenvalue counts do not depend on the vertex order)."""
    n = k1 + k2 + p + 2
    c1, c2 = k1, n - 1
    edges = [(i, c1) for i in range(k1)]
    chain = [c1] + list(range(k1 + 1, k1 + 1 + p)) + [c2]
    edges += [(chain[i], chain[i + 1]) for i in range(len(chain) - 1)]
    edges += [(i, c2) for i in range(k1 + 1 + p, k1 + 1 + p + k2)]
    return n, edges


def degrees(n, edges):
    d = [0] * n
    for u, v in edges:
        d[u] += 1
        d[v] += 1
    return d


def check_tree(n, edges, expected_degrees):
    """Sanity: simple graph, n-1 edges, connected (so a tree), expected degree multiset."""
    require(len(edges) == n - 1, "wrong edge count")
    require(all(0 <= u < n and 0 <= v < n for u, v in edges), "vertex out of range")
    require(len({frozenset(e) for e in edges}) == n - 1 and all(u != v for u, v in edges), "graph is not simple")
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    seen, stack = {0}, [0]
    while stack:
        for w in adj[stack.pop()]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    require(len(seen) == n, "graph is not connected")
    require(sorted(degrees(n, edges)) == sorted(expected_degrees), "unexpected degrees")


def dc_expected_degrees(k1, k2, p):
    return [k1 + 1, k2 + 1] + [2] * p + [1] * (k1 + k2)


def laplacian(n, edges):
    L = [[0] * n for _ in range(n)]
    for u, v in edges:
        L[u][v] -= 1
        L[v][u] -= 1
        L[u][u] += 1
        L[v][v] += 1
    return L


# ======================================================================================
# Randic index: rational enclosure
# ======================================================================================
def inv_sqrt_enclosure(m):
    """Rationals lo <= 1/sqrt(m) <= hi (equality iff m is a perfect square)."""
    s = isqrt(m)
    if s * s == m:
        return Fr(1, s), Fr(1, s)
    S = isqrt(m * 10 ** (2 * R_DIGITS))      # S < 10^D sqrt(m) < S + 1  (sqrt(m) irrational)
    return Fr(10 ** R_DIGITS, S + 1), Fr(10 ** R_DIGITS, S)


def randic_enclosure(n, edges):
    d = degrees(n, edges)
    lo = hi = Fr(0)
    for u, v in edges:
        l, h = inv_sqrt_enclosure(d[u] * d[v])
        lo += l
        hi += h
    return lo, hi


def randic_float_closed_form_dc(k1, k2, p):
    """Float value of R(DC(k1,k2,p)) from its closed form (sanity check of the constructor)."""
    r = k1 / sqrt(k1 + 1) + k2 / sqrt(k2 + 1)
    if p == 0:
        return r + 1 / sqrt((k1 + 1) * (k2 + 1))
    return r + 1 / sqrt(2 * (k1 + 1)) + 1 / sqrt(2 * (k2 + 1)) + (p - 1) / 2


# ======================================================================================
# Method A: Sylvester inertia of L - tI via exact LDL^T
# ======================================================================================
def count_eigs_below(L, t):
    """Number of eigenvalues of the symmetric integer matrix L that are < t, i.e. the number
    of negative pivots of the LDL^T factorisation of L - tI (Sylvester's law of inertia).
    Exact when t is a Fraction (a zero pivot raises ZeroDivisionError).  With a float t the
    same code is only a heuristic, used for locating a."""
    exact = isinstance(t, Fr)
    n = len(L)
    A = [[L[i][j] - t if i == j else L[i][j] for j in range(n)] for i in range(n)]
    neg = 0
    for k in range(n):
        piv = A[k][k]
        if piv == 0:
            if exact:
                raise ZeroDivisionError("zero pivot")
            piv = 1e-300
        if piv < 0:
            neg += 1
        rowk = A[k]
        for i in range(k + 1, n):
            aik = A[i][k]
            if aik != 0:
                f = aik / piv
                Ai = A[i]
                for j in range(k + 1, n):
                    if rowk[j] != 0:
                        Ai[j] -= f * rowk[j]
    return neg


def a_float_guess(L):
    """Float guess for a by bisection on the (float) eigenvalue count.  Never trusted."""
    lo, hi = 0.0, float(len(L))
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            break
        if count_eigs_below(L, mid) >= 2:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def a_by_inertia(L):
    """Certified lo <= a < hi via two exact inertia counts.  Returns (lo, hi, (c_lo, c_hi))."""
    g = Fr(a_float_guess(L))
    for D in (A_DIGITS, A_DIGITS - 3, A_DIGITS - 6, A_DIGITS - 9):
        m = (g.numerator * 10 ** D) // g.denominator          # floor(g * 10^D)
        lo, hi = Fr(m - 1, 10 ** D), Fr(m + 2, 10 ** D)
        if lo <= 0:
            continue
        try:
            c_lo, c_hi = count_eigs_below(L, lo), count_eigs_below(L, hi)
        except ZeroDivisionError:
            continue
        if c_lo == 1 and c_hi >= 2:
            return lo, hi, (c_lo, c_hi)
    raise RuntimeError("Method A could not certify a")


# ======================================================================================
# Method B: exact characteristic polynomial + Sturm sequence
# ======================================================================================
def matmul(X, Y):
    n = len(X)
    Yt = list(zip(*Y))
    return [[sum(a * b for a, b in zip(X[i], Yt[j])) for j in range(n)] for i in range(n)]


def charpoly(L):
    """Integer coefficients c[0..n] of det(xI - L) = sum c[i] x^i (Faddeev-LeVerrier over Z):
    M_0 = 0, M_k = L M_{k-1} + c_{n-k+1} I, c_{n-k} = -tr(L M_k)/k (an exact division)."""
    n = len(L)
    c = [0] * (n + 1)
    c[n] = 1
    M = [[0] * n for _ in range(n)]
    for k in range(1, n + 1):
        M = matmul(L, M)
        for i in range(n):
            M[i][i] += c[n - k + 1]
        tr = sum(L[i][j] * M[j][i] for i in range(n) for j in range(n))
        if tr % k:
            raise ArithmeticError("Faddeev-LeVerrier: inexact division")
        c[n - k] = -(tr // k)
    LM = matmul(L, M)                                   # Cayley-Hamilton: L M_n + c_0 I = 0
    require(all(LM[i][j] + (c[0] if i == j else 0) == 0 for i in range(n) for j in range(n)), "Cayley-Hamilton check failed")
    return c


def poly_trim(p):
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def poly_primitive(p):
    """The positive rational multiple of p with coprime integer coefficients."""
    p = [Fr(x) for x in p]
    den = 1
    for x in p:
        den = den * x.denominator // gcd(den, x.denominator)
    ints = [int(x * den) for x in p]
    g = 0
    for v in ints:
        g = gcd(g, v)
    return poly_trim([v // g for v in ints] if g else ints)


def poly_rem(a, b):
    """Remainder of a modulo b over Q (coefficient lists, lowest degree first)."""
    r = [Fr(x) for x in a]
    db, lb = len(b) - 1, Fr(b[-1])
    while len(r) - 1 >= db and any(r):
        f = r[-1] / lb
        s = len(r) - 1 - db
        for i in range(db + 1):
            r[s + i] -= f * b[i]
        r.pop()
        poly_trim(r)
    return r if r else [Fr(0)]


def sturm_chain(q):
    """Sturm sequence s0 = q, s1 = q', s_{i+1} = -rem(s_{i-1}, s_i), each scaled by a positive
    constant to a primitive integer polynomial (positive scaling preserves the Sturm property)."""
    S = [poly_primitive(q), poly_primitive([i * q[i] for i in range(1, len(q))])]
    while True:
        r = poly_rem(S[-2], S[-1])
        if not any(r):
            return S
        S.append(poly_primitive([-x for x in r]))


def sign_at(p, t):
    """Exact sign of the integer polynomial p at the rational t (homogenised Horner)."""
    num, den = t.numerator, t.denominator
    v, pw = p[-1], 1
    for c in reversed(p[:-1]):
        pw *= den
        v = v * num + c * pw
    return (v > 0) - (v < 0)


def variations(signs):
    s = [x for x in signs if x != 0]
    return sum(1 for i in range(len(s) - 1) if s[i] != s[i + 1])


def V(S, t):
    return variations([sign_at(s, t) for s in S])


def V_minus_inf(S):
    return variations([(1 if s[-1] > 0 else -1) * (-1) ** (len(s) - 1) for s in S])


def a_by_sturm(L):
    """Certified lo < a < hi with no floating point.  Returns (lo, hi, info, charpoly)."""
    n = len(L)
    c = charpoly(L)
    require(c[0] == 0, "det(L) is not zero")
    q = c[1:]                              # q(x) = det(xI - L) / x
    if q[0] == 0:
        raise RuntimeError("0 is a multiple eigenvalue: graph disconnected")
    S = sturm_chain(q)
    vinf = V_minus_inf(S)
    lo, hi = Fr(0), Fr(2 * n + 1, 2)       # all Laplacian eigenvalues lie in [0, n]
    vlo, vhi = V(S, lo), V(S, hi)
    require(vinf - vlo == 0 and vlo - vhi >= 1, "invalid initial Sturm bracket")
    width = Fr(1, 10 ** (B_DIGITS + 1))
    while hi - lo > width:                 # invariant: no root <= lo, some root in (lo, hi)
        mid = (lo + hi) / 2
        while sign_at(q, mid) == 0:        # only possible at an integer mid
            mid = (mid + hi) / 2
        vm = V(S, mid)
        if vlo - vm >= 1:
            hi = mid
        else:
            lo, vlo = mid, vm
    # round outward to the decimal grid and re-certify from scratch
    D = 10 ** B_DIGITS
    lo_d = Fr((lo.numerator * D) // lo.denominator, D)
    hi_d = Fr(-((-hi.numerator * D) // hi.denominator), D)
    while sign_at(q, lo_d) == 0:
        lo_d -= Fr(1, D)
    while sign_at(q, hi_d) == 0:
        hi_d += Fr(1, D)
    roots_below = vinf - V(S, lo_d)        # distinct roots of q in (-inf, lo_d]
    roots_inside = V(S, lo_d) - V(S, hi_d)  # distinct roots of q in (lo_d, hi_d)
    if not (roots_below == 0 and roots_inside >= 1):
        raise RuntimeError("Method B could not certify a")
    return lo_d, hi_d, (roots_below, roots_inside, len(S)), c


def charpoly_sympy(L):
    x = sympy.Symbol("x")
    P = sympy.Matrix(L).charpoly(x)
    return [int(v) for v in reversed(P.all_coeffs())]


# ======================================================================================
# Method C (paths only): a(P_n) = 2 - 2 cos(pi/n) with rational bounds
# ======================================================================================
def alt_series_enclosure(terms):
    """For an alternating series with terms decreasing in absolute value, the sum lies
    between any two consecutive partial sums."""
    s, partial = Fr(0), []
    for t in terms:
        s += t
        partial.append(s)
    return min(partial[-2:]), max(partial[-2:])


def atan_inv(m, N):          # arctan(1/m), m >= 2
    x = Fr(1, m)
    return alt_series_enclosure([Fr((-1) ** k, 2 * k + 1) * x ** (2 * k + 1) for k in range(N)])


def pi_enclosure():
    a5, a239 = atan_inv(5, 40), atan_inv(239, 15)          # Machin: pi = 16 atan(1/5) - 4 atan(1/239)
    lo, hi = 16 * a5[0] - 4 * a239[1], 16 * a5[1] - 4 * a239[0]
    D = 10 ** 50                                           # round outward to keep numbers small
    return Fr((lo.numerator * D) // lo.denominator, D), Fr(-((-hi.numerator * D) // hi.denominator), D)


def cos_enclosure(x, N=30):  # 0 < x < 1
    terms, t = [], Fr(1)
    for k in range(N):
        if k:
            t = -t * x * x / ((2 * k - 1) * (2 * k))
        terms.append(t)
    return alt_series_enclosure(terms)


def a_path_closed_form(n):
    plo, phi = pi_enclosure()
    xlo, xhi = plo / n, phi / n                            # pi/n in [xlo, xhi], inside (0, 1)
    require(0 < xlo < xhi < 1, "cosine series argument outside certified domain")
    cos_lo = cos_enclosure(xhi)[0]                         # cos is decreasing on [0, pi]
    cos_hi = cos_enclosure(xlo)[1]
    return 2 - 2 * cos_hi, 2 - 2 * cos_lo


# ======================================================================================
# Certification of one graph
# ======================================================================================
def certify(n, edges, use_sympy, is_path=False, dc=None):
    if is_path:
        check_tree(n, edges, [1, 1] + [2] * (n - 2))
    else:
        check_tree(n, edges, dc_expected_degrees(*dc))
    L = laplacian(n, edges)
    Rlo, Rhi = randic_enclosure(n, edges)
    Rf = (n - 3) / 2 + sqrt(2) if is_path else randic_float_closed_form_dc(*dc)
    require(Rlo - Fr(1, 10 ** 9) <= Fr(Rf) <= Rhi + Fr(1, 10 ** 9), "R closed-form sanity check failed")
    aA_lo, aA_hi, certA = a_by_inertia(L)
    aB_lo, aB_hi, certB, cp = a_by_sturm(L)
    consistent = max(aA_lo, aB_lo) <= min(aA_hi, aB_hi)    # both enclose the same number
    rec = dict(n=n, R=(Rlo, Rhi), aA=(aA_lo, aA_hi), aB=(aB_lo, aB_hi), certA=certA, certB=certB,
               RaA=(Rlo * aA_lo, Rhi * aA_hi), RaB=(Rlo * aB_lo, Rhi * aB_hi))
    if use_sympy:
        rec["sympy_charpoly_equal"] = (charpoly_sympy(L) == cp)
        consistent = consistent and rec["sympy_charpoly_equal"]
    if is_path:
        aC = a_path_closed_form(n)
        rec["aC"] = aC
        consistent = consistent and max(aA_lo, aB_lo, aC[0]) <= min(aA_hi, aB_hi, aC[1])
    rec["consistent"] = consistent
    return rec


def double_comets(n):
    """All (k1, k2, p) with k1, k2 >= 2, p >= 0, k1 + k2 + p + 2 = n (ordered pairs)."""
    return [(k1, n - 2 - p - k1, p) for p in range(0, n - 5) for k1 in range(2, n - 3 - p)]


# ======================================================================================
# Output helpers (directed decimal rounding of exact rationals)
# ======================================================================================
def dec(x, d=SHOW, up=False, signed=False):
    """Decimal string of the rational x with d decimals, rounded down (up=False) or up."""
    x = Fr(x)
    m = -((-x.numerator * 10 ** d) // x.denominator) if up else (x.numerator * 10 ** d) // x.denominator
    s = str(abs(m)).rjust(d + 1, "0")
    return ("-" if m < 0 else "+" if signed else "") + s[:-d] + "." + s[-d:]


def ivl(I, d=SHOW, signed=False):
    return "[" + dec(I[0], d, signed=signed) + ", " + dec(I[1], d, up=True, signed=signed) + "]"


def sgn(x, d=SHOW):          # signed, rounded toward -inf (used for lower bounds)
    return dec(x, d, signed=True)


# ======================================================================================
# Driver
# ======================================================================================
def run_n(n, use_sympy, out, mode):
    """mode 'theorem': every DC must be certified above P_n.
       mode 'table'  : classify each DC as ABOVE / BELOW / undecided (both methods)."""
    P = certify(*path_graph(n), use_sympy, is_path=True)
    out(f"n = {n}")
    out(f"  P_{n}:  R in {ivl(P['R'])}")
    out(f"         a in {ivl(P['aA'])} (A: inertia counts {P['certA']})")
    out(f"         a in {ivl(P['aB'])} (B: Sturm, roots below/inside = {P['certB'][:2]})")
    out(f"         a in {ivl(P['aC'])} (C: 2-2cos(pi/n), rational pi & cos bounds)")
    out(f"         R*a(P_{n}) in {ivl(P['RaA'])} (A)   {ivl(P['RaB'])} (B)")
    rows, ordered_ok, ordered_total, consistent = [], 0, 0, P["consistent"]
    recs = {}
    for (k1, k2, p) in double_comets(n):
        n2, e = double_comet(k1, k2, p)
        require(n2 == n, "double comet has wrong order")
        recs[(k1, k2, p)] = certify(n, e, use_sympy, dc=(k1, k2, p))
    for (k1, k2, p), r in recs.items():
        consistent = consistent and r["consistent"]
        gapA = r["RaA"][0] - P["RaA"][1]     # lower(DC) - upper(P), Method A
        gapB = r["RaB"][0] - P["RaB"][1]     # lower(DC) - upper(P), Method B
        revA = P["RaA"][0] - r["RaA"][1]     # lower(P) - upper(DC), Method A
        revB = P["RaB"][0] - r["RaB"][1]
        r.update(gapA=gapA, gapB=gapB, revA=revA, revB=revB,
                 above=(gapA > 0 and gapB > 0), below=(revA > 0 and revB > 0))
        ordered_total += 1
        ordered_ok += r["above"]
    # one row per isomorphism class (DC(k1,k2,p) = DC(k2,k1,p)); both orders were certified
    for (k1, k2, p), r in recs.items():
        if k1 > k2:
            continue
        twin = recs[(k2, k1, p)]
        rows.append((k1, k2, p, r, twin))
    rows.sort(key=lambda x: x[3]["RaA"][0])
    if mode == "theorem":
        out(f"  {'DC(k1,k2,p)':<13}{'R*a(DC) >= (A)':<18}{'R*a(DC) >= (B)':<18}"
            f"{'gap = R*a(DC) - R*a(P_n) >=':<30}{'DC > P_n certified'}")
    else:
        out(f"  {'DC(k1,k2,p)':<13}{'R*a(DC) in (A)':<36}{'R*a(DC) - R*a(P_n) in (A)':<39}{'status (A and B)'}")
    for k1, k2, p, r, twin in rows:
        tag = f"DC({k1},{k2},{p})"
        if mode == "theorem":
            gap = min(r["gapA"], r["gapB"], twin["gapA"], twin["gapB"])
            ok = r["above"] and twin["above"]
            out(f"  {tag:<13}{dec(r['RaA'][0]):<18}{dec(r['RaB'][0]):<18}{sgn(gap):<30}{ok}")
        else:
            diff = (r["RaA"][0] - P["RaA"][1], r["RaA"][1] - P["RaA"][0])   # encloses R*a(DC) - R*a(P_n)
            status = ("DC > P_n" if (r["above"] and twin["above"]) else
                      "DC < P_n  (reversal)" if (r["below"] and twin["below"]) else "undecided")
            out(f"  {tag:<13}{ivl(r['RaA']):<36}{ivl(diff, signed=True):<39}{status}")
    best = rows[0]
    # Sorting lower endpoints alone is not a certificate of the minimum.
    require(all(best[3][key][1] < row[3][key][0]
                for row in rows[1:] for key in ("RaA", "RaB")),
            "minimum among proper double comets is not isolated")
    out(f"  -> {len(rows)} isomorphism classes, {ordered_total} ordered triples (k1,k2,p); "
        f"certified DC > P_n for {ordered_ok}/{ordered_total}; methods A/B/C consistent: {consistent}"
        + (f"; sympy charpoly agrees: {all(x.get('sympy_charpoly_equal', True) for x in list(recs.values()) + [P])}"
           if use_sympy else ""))
    out(f"  -> smallest R*a among double comets: DC({best[0]},{best[1]},{best[2]}) "
        f"in {ivl(best[3]['RaA'])}")
    out("")
    return P, recs, ordered_ok == ordered_total, consistent


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--n", type=int, nargs="+", help="certify comparison tables for these n (n >= 6)")
    ap.add_argument("--no-sympy", action="store_true", help="skip the optional sympy cross-check")
    args = ap.parse_args()
    if args.n and any(n < 6 for n in args.n):
        ap.error("--n values must be at least 6")
    use_sympy = (sympy is not None) and not args.no_sympy
    out = lambda s="": print(s, flush=True)
    t0 = time.time()
    out("Certified verification: R*a of double comets vs. the path P_n")
    out("  R = Randic index, a = algebraic connectivity; DC(k1,k2,p): centres with k1, k2 >= 2")
    out("  pendant vertices joined by a path with p >= 0 internal vertices (n = k1+k2+p+2).")
    out("  Intervals are exact rational enclosures, printed with lower ends rounded down and")
    out("  upper ends rounded up.  A = Sylvester inertia (exact LDL^T of L - tI); B = Sturm")
    out("  sequence of det(xI - L)/x over Z (no floating point); C = 2-2cos(pi/n) (paths only).")
    out(f"  python {platform.python_version()}; sympy cross-check: "
        + (f"on (sympy {sympy.__version__})" if use_sympy else "off"))
    out("")
    if args.n:
        good = True
        for n in args.n:
            _, recs, _, cons = run_n(n, use_sympy, out, "table")
            good = good and cons and all(r["above"] or r["below"] for r in recs.values())
        out(f"runtime: {time.time() - t0:.1f} s")
        sys.exit(0 if good else 1)

    out("=== THEOREM: for n = 10, 11, 12, 13 every proper double comet has R*a > R*a(P_n) ===")
    out("")
    thm = True
    for n in (10, 11, 12, 13):
        _, _, ok, cons = run_n(n, use_sympy, out, "theorem")
        thm = thm and ok and cons
    out(f"THEOREM (n = 10..13) CERTIFIED: {thm}")
    out("")
    out("=== n = 14: the inequality reverses for DC(2,2,8) ===")
    out("")
    P, recs, _, cons14 = run_n(14, use_sympy, out, "table")
    r = recs[(2, 2, 8)]
    rev = r["below"] and cons14 and all(
        x["above"] for key, x in recs.items() if key != (2, 2, 8))
    out(f"  R*a(DC(2,2,8)) <= {dec(r['RaA'][1], up=True)} (A), {dec(r['RaB'][1], up=True)} (B)"
        f"  <  R*a(P_14) >= {dec(P['RaA'][0])} (A), {dec(P['RaB'][0])} (B)")
    out(f"  certified gap R*a(P_14) - R*a(DC(2,2,8)) >= {sgn(min(r['revA'], r['revB']))}")
    out(f"REVERSAL AT n = 14 FOR DC(2,2,8) CERTIFIED: {rev}")
    out("")
    out(f"runtime: {time.time() - t0:.1f} s")
    out(f"ALL CERTIFIED: {thm and rev}")
    sys.exit(0 if (thm and rev) else 1)


if __name__ == "__main__":
    main()
