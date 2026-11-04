# Exact verification

## Primary interval verifier

`certificate/verify_theorem.py` constructs all ordered triples `(k1,k2,p)` of the specified orders with `k1,k2 >= 2`, `p >= 0`, and `k1+k2+p+2=n`. It checks graph simplicity, connectedness, order and degree multisets. Swapping the two star sizes produces the same isomorphism class; both orders are verified.

For each edge, integer square roots give rational lower and upper bounds on `1/sqrt(d(u)d(v))`. Summing gives an enclosure for `R`.

Two methods enclose the algebraic connectivity of **each graph**:

1. Exact rational elimination of `L-tI`, with Sylvester inertia counts at rational endpoints. Floating point proposes endpoints but cannot certify them. Zero pivots abort that attempt. Non-integer rational endpoints cannot be roots of the monic integer leading-principal-minor polynomials.
2. Exact Faddeev–LeVerrier characteristic polynomial, a Cayley–Hamilton check, and a Sturm sequence to isolate the smallest root after dividing out the simple zero eigenvalue. This method uses no floating point.

The third spectral check applies **only to paths**: the formula `a(P_n)=2-2cos(pi/n)` is enclosed using rational Machin-series bounds for pi and alternating cosine-series bounds. When SymPy is installed, its characteristic polynomial is compared with the first implementation's polynomial. This is a polynomial cross-check, not a third eigenvalue method for every graph.

All inequalities in the theorem are decided from disjoint rational product intervals. The printed decimals are rounded outward. The minimum within each proper-DC family is checked by comparing the candidate's upper endpoint with every other class's lower endpoint, under both methods; sorting lower endpoints alone is not used as a proof. The default run also checks that at order 14 the remaining proper-DC classes are above the path.

## Independent determinant verifier

`certificate/verify_principal_minors.py` imports no code from the primary verifier. It uses a different labelling, a fraction-free Bareiss determinant with exact-division checks, and the separating rational constants in the README.

Let `J` be the all-ones matrix. For a connected graph on `n` vertices the spectrum of `L+J` is `n, lambda_2,...,lambda_n`. For rational `0<t<n`:

- Positive leading principal minors of `L+J-tI` imply `a>t`, by Sylvester's criterion.
- A negative determinant implies a negative eigenvalue of `L+J-tI`. Since `n-t>0`, this gives `a<t`.

For a rational target `c`, take `t=c/R_lower` to prove `R*a>c`, or `t=c/R_upper` to prove `R*a<c`. Denominators are cleared by a **positive** integer before determinant calculation, preserving signs. The script verifies every radical enclosure by squaring rational endpoints. It checks 57 proper-DC classes at orders 10–13, their four path comparisons, and two comparisons at order 14.

This certificate uses no floating point, no numerical eigenvalue solver, no characteristic polynomial, and no imports from the primary verifier. Independence here means a separate implementation and mathematical criterion, not an external human review.

## Commands

```sh
python3 -I certificate/verify_theorem.py --no-sympy
python3 -I -O certificate/verify_theorem.py --no-sympy
python3 -I -O certificate/verify_principal_minors.py
```

Optional polynomial cross-check:

```sh
python3 -m pip install -r requirements-audit.txt
python3 -I certificate/verify_theorem.py
```

Every mandatory check uses explicit exceptions, so optimization cannot remove it. CI runs the exact verifiers in normal and optimized modes. This is a computer-assisted exact certificate, **not a Lean formalization**.
