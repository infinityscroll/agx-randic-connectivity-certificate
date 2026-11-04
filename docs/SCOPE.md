# Statement and interpretation

The unconditional theorem in this repository compares the path with every `DC(k1,k2,p)` satisfying `k1,k2 >= 2`, for each `n=10,11,12,13`. It also proves one reversed comparison at `n=14`. The parameter `p` counts internal vertices, not edges, on the path between the centres.

Both the [2007 original](https://match.pmf.kg.ac.rs/electronic_versions/Match58/n1/match58n1_83-102.pdf) (Conjecture 4, p. 98) and the [2012 restatement](https://match.pmf.kg.ac.rs/electronic_versions/Match68/n3/match68n3_843-849.pdf) (Conjecture 1.2, p. 845) propose paths through order 9 and balanced double comets from order 10. Their parenthetical descriptions specify two equal-sized stars joined by a path, without an explicit minimum star size.

Under the proper convention, both centres have degree at least three. The theorem excludes every proposed balanced minimizer for orders 10–13 by exhibiting a connected graph with a smaller value. No proof of global path optimality is needed for this conditional refutation.

Under the inclusive convention, `DC(1,1,n-4)` is isomorphic to the path `P_n`. Our comparisons do not contradict an assertion that *some* member of that inclusive family minimizes the product. A partition of a statement into two ranges does not logically require the named families to be disjoint.

A later paper by Aouchiche, Hansen, Hertz and Sellal uses a balanced-double-comet convention with diameter at most `n-3` ([GERAD G-2015-31](https://www.gerad.ca/fr/papers/G-2015-31.pdf), Section 2). That is evidence for a proper convention in related literature, but it does not remove the ambiguity in the 2007 wording.

Accordingly this release does **not** use “Conjecture 4 refuted as stated” as its headline. It establishes a finite comparison theorem and a conditional refutation of the product clause. It does not settle the full asymptotic minimization problem, the sum clause, or the quotient clause.
