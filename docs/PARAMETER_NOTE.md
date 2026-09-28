# Thesis-parameter note

The Appendix-B program is a general parameterized simulation program, not a single fixed experiment. For reproducing a thesis result, the input parameters must be set to the values used in that result.

The frozen thesis-validation preset therefore uses `collision_alpha=3`, following thesis Sec. 2.3. It also uses first-order weighting, following Sec. 3.1. The portfolio validation wrapper defaults to three repetitions and increases to five only if stochastic variability makes a representative result ambiguous.
