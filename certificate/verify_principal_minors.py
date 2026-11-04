#!/usr/bin/env python3
"""Independent exact R*a comparisons; no imports from the original certificate.

For L a connected graph Laplacian, L+J has eigenvalues n, lambda_2,...,lambda_n.
Sylvester's criterion certifies lambda_2 > t from positive leading principal
minors of L+J-tI (with t<n). A negative determinant certifies lambda_2 < t.
All arithmetic is integer/rational. No numerical eigenvalues are computed.
"""
from fractions import Fraction as Q
from math import isqrt
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def det(a):
    """Fraction-free Bareiss determinant, with row swaps and checked divisions."""
    a = [list(row) for row in a]
    n = len(a)
    prev, sign = 1, 1
    for k in range(n - 1):
        if not a[k][k]:
            rows = [i for i in range(k + 1, n) if a[i][k]]
            if not rows:
                return 0
            j = rows[0]
            a[k], a[j] = a[j], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                value = pivot * a[i][j] - a[i][k] * a[k][j]
                need(value % prev == 0, 'Bareiss division not exact')
                a[i][j] = value // prev
            a[i][k] = 0
        prev = pivot
    return sign * a[-1][-1]


def graph(n, leaves=None):
    """Build path or DC using centers at positions 0 and p+1 on the spine."""
    if leaves is None:
        edges = [(i, i + 1) for i in range(n - 1)]
    else:
        left, right = leaves
        p = n - left - right - 2
        need(left >= 2 and right >= 2 and p >= 0, 'Bad double-comet parameters')
        edges = [(i, i + 1) for i in range(p + 1)]
        edges += [(0, p + 2 + i) for i in range(left)]
        edges += [(p + 1, p + 2 + left + i) for i in range(right)]
    need(len(edges) == n - 1, 'Wrong edge count')
    degree = [0] * n
    L = [[0] * n for _ in range(n)]
    adjacent = [set() for _ in range(n)]
    for i, j in edges:
        need(0 <= i < n and 0 <= j < n and i != j and j not in adjacent[i], 'Bad edge')
        adjacent[i].add(j)
        adjacent[j].add(i)
        degree[i] += 1
        degree[j] += 1
        L[i][j] -= 1
        L[j][i] -= 1
        L[i][i] += 1
        L[j][j] += 1
    seen, todo = {0}, [0]
    while todo:
        u = todo.pop()
        for v in adjacent[u] - seen:
            seen.add(v)
            todo.append(v)
    need(len(seen) == n, 'Disconnected graph')
    if leaves is not None:
        p = n - sum(leaves) - 2
        need(sorted(degree) == sorted([leaves[0]+1, leaves[1]+1]+[1]*sum(leaves)+[2]*p), 'Bad DC degrees')
    else:
        need(sorted(degree) == sorted([1, 1] + [2]*(n-2)), 'Bad path degrees')
    scale = 10**12
    Rlo, Rhi = Q(0), Q(0)
    for u, v in edges:
        m = degree[u] * degree[v]
        floor = isqrt(m * scale**2)
        lo, hi = Q(scale, floor+1), Q(scale, floor)
        need(0 < lo <= hi and m*lo*lo <= 1 <= m*hi*hi, 'Bad radical bound')
        Rlo += lo
        Rhi += hi
    return L, Rlo, Rhi


def comparison(n, leaves, target, above):
    L, Rlo, Rhi = graph(n, leaves)
    t = target / (Rlo if above else Rhi)
    need(0 < t < n, 'Invalid spectral threshold')
    a, b = t.numerator, t.denominator
    M = [[b*(L[i][j]+1) - (a if i == j else 0) for j in range(n)] for i in range(n)]
    if above:
        signs = [det([row[:k] for row in M[:k]]) > 0 for k in range(1, n+1)]
        need(all(signs), f'Could not prove R*a > {target}: n={n}, leaves={leaves}')
    else:
        need(det(M) < 0, f'Could not prove R*a < {target}: n={n}, leaves={leaves}')
    return {'n':n, 'leaves':leaves, 'relation':'>' if above else '<', 'target':str(target),
            'spectral_threshold':str(t), 'method':'all leading minors positive' if above else 'determinant negative'}


def main():
    # Elementary determinant sanity checks include zero, row swap, and nonunit pivots.
    need(det([[0,1],[1,0]]) == -1, 'Swap regression')
    need(det([[1,2],[2,4]]) == 0, 'Singular regression')
    need(det([[2,1,3],[1,4,2],[3,2,5]]) == 3, 'Bareiss regression')
    targets = {10:Q(49,100), 11:Q(44,100), 12:Q(405,1000), 13:Q(373,1000)}
    records, counts = [], {}
    for n, target in targets.items():
        records.append(comparison(n, None, target, False))
        count = 0
        for left in range(2, n):
            for right in range(left, n):
                if left + right + 2 <= n:
                    records.append(comparison(n, (left,right), target, True))
                    count += 1
        counts[n] = count
    need(counts == {10:9,11:12,12:16,13:20}, 'Enumeration count mismatch')
    records.append(comparison(14, None, Q(346,1000), True))
    records.append(comparison(14, (2,2), Q(346,1000), False))
    print(json.dumps({'all_certified':True,'classes':counts,'comparisons':records},indent=2))


if __name__ == '__main__':
    main()
