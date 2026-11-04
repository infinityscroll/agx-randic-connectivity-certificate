# Paths versus proper double comets for Randić index times algebraic connectivity

[![Exact certificates](https://github.com/infinityscroll/agx-randic-connectivity-certificate/actions/workflows/verify.yml/badge.svg)](https://github.com/infinityscroll/agx-randic-connectivity-certificate/actions/workflows/verify.yml)

**Jitendra Prajapati · Independent · [ORCID](https://orcid.org/0009-0008-7493-0311)**

Exact arithmetic proves that, for each order **10, 11, 12, and 13**, the path has a smaller product of Randić index and algebraic connectivity than **every proper double comet** of that order, balanced or unbalanced. At order **14**, the balanced double comet with two leaves at each end beats the path.

This bears on the product clause of Aouchiche–Hansen–Zheng (2007), Conjecture 4, restated by Li–Shi–Wang (2012), Conjecture 1.2. **It refutes the clause if “double comet” excludes paths.** The printed parenthetical definition does not specify a minimum star size; if one-leaf stars are allowed, a path is itself a balanced double comet, and this result does not refute that interpretation. The original path/double-comet case split alone does not resolve this ambiguity. See [scope](docs/SCOPE.md) and [literature review](docs/PRIORITY.md).

## Read and reproduce

- [Mathematical note (PDF)](paper/note.pdf) and [LaTeX source](paper/note.tex).
- [Proof and verifier details](docs/VERIFICATION.md).
- [Exact interval output](results/exact_certificate.txt) and [independent determinant certificate](results/principal_minors.json).
- [Numerical observations](docs/NUMERICAL_EVIDENCE.md), separated from the theorem.

Python 3.10 or later; the two exact verifiers use only the standard library:

```sh
python3 -I certificate/verify_theorem.py --no-sympy
python3 -I certificate/verify_principal_minors.py
```

The first command must end with `ALL CERTIFIED: True`. The second must report `"all_certified": true`. Both exit nonzero on a failed certificate. Checks remain active with `python3 -I -O`.

## Precise finite theorem

Let `DC(k1,k2,p)` consist of two centres connected by a path with `p` internal vertices and with `k1` and `k2` pendant leaves at the centres. Here **`k1,k2 >= 2`**, `p >= 0`, and `n = k1+k2+p+2`. We call these *proper double comets*.

Write `R(G) = sum_{uv in E(G)} 1/sqrt(d(u)d(v))` and let `a(G)` be the second-smallest Laplacian eigenvalue. The independent verifier proves the following strict separating bounds without computing numerical eigenvalues:

| Order | Proper DC classes checked | Exact separation, for every such DC |
|---:|---:|---|
| 10 | 9 | `R(P10)a(P10) < 49/100 < R(DC)a(DC)` |
| 11 | 12 | `R(P11)a(P11) < 11/25 < R(DC)a(DC)` |
| 12 | 16 | `R(P12)a(P12) < 81/200 < R(DC)a(DC)` |
| 13 | 20 | `R(P13)a(P13) < 373/1000 < R(DC)a(DC)` |

It also proves `R(DC(2,2,8))a(DC(2,2,8)) < 173/500 < R(P14)a(P14)`.

The interval verifier separately certifies that `DC(2,2,n-6)` is the unique minimum **within the proper-double-comet family** for each `n = 10,...,14`, up to swapping the two centres. These are finite statements, not an assertion for all larger orders.

## Limits of the result

- This theorem does not establish the global minimizer over all connected graphs, or over all trees.
- Saved floating-point scans report the path as the minimum among 11,716,571 connected graphs at order 10, and a path/double-comet pattern among trees through order 22. Those are numerical observations, not exact global certificates. Only retained candidates were replayed during this review.
- The incomplete connected-graph scan at order 11 is excluded.
- No earlier result with this precise finite product comparison was located in the reviewed sources as of November 5, 2026. This is a search result, **not a priority guarantee or proof that the problem remained open**.
- This addresses one product clause, not the other clauses of Conjecture 4, and not a general classification of extremal graphs.

## Sources and credit

1. M. Aouchiche, P. Hansen, M. Zheng, *Variable Neighborhood Search for Extremal Graphs. 19. Further Conjectures and Results about the Randić Index*, MATCH 58 (2007), 83–102, Conjecture 4, p. 98. [Publisher PDF](https://match.pmf.kg.ac.rs/electronic_versions/Match58/n1/match58n1_83-102.pdf).
2. X. Li, Y. Shi, L. Wang, *On a Relation Between Randić Index and Algebraic Connectivity*, MATCH 68 (2012), 843–849, Conjecture 1.2, p. 845. [Publisher PDF](https://match.pmf.kg.ac.rs/electronic_versions/Match68/n3/match68n3_843-849.pdf).
3. J. J. Gebendorfer, *An infinite tadpole counterexample family to the Randić index–algebraic connectivity conjecture* (2026). This concerns the **quotient R/a**, a different clause. [DOI](https://doi.org/10.5281/zenodo.21827918).

See [provenance](docs/PROVENANCE.md) for computational assistance and review history. Code is MIT licensed; the note and documentation are CC BY 4.0. Third-party papers are linked, not redistributed. A GitHub release is a reproducible research artifact, not journal publication or human peer review.
