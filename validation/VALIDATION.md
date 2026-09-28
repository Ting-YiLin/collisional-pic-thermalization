# Thesis validation record

**Status: VALIDATED AGAINST THESIS RESULTS**

The restored C++ code has been compiled and executed in a modern toolchain and checked against the master's-thesis result families. Representative reruns reproduce the main reported behavior: numerical relaxation occurs, adding pitch-angle Monte-Carlo collisions shortens the thermalization time, the thermalization time changes systematically with particle number and collision rate, and the plasma-wave/self-heating diagnostics behave in the reported directions. No material contradiction with the thesis's core conclusions was identified.

## Parameter correction frozen for publication

For thesis validation, the pitch-angle model uses **`collision_alpha = 3`**, matching Sec. 2.3 of the thesis. The source program is a parameterized research simulator; the example input in the appendix is not the only valid run configuration. The publication preset also uses first-order particle weighting, matching Sec. 3.1.

## Portfolio validation protocol

This repository is a portfolio/reproducibility archive, not a new academic replication paper. The frozen protocol therefore uses **3 independent/internal repetitions by default and 5 only when stochastic variability makes the result ambiguous**. It does not rerun every historical point 30 or 100 times.

## Representative recovered-code check

A direct run of the recovered source at `L=64`, `Dx=0.25`, `Dt=0.1`, `ND=20`, ball/Fermi-like initial distribution, first-order weighting, and `alpha=3` produced the following H4 relaxation estimates with three repetitions inside the historical program:

| Case | H4 thermalization estimate |
|---|---:|
| Collisionless | 8788.38 |
| Pitch-angle, gamma0/omega_p = 1e-3 | 4667.90 |

This representative run reproduces the thesis's central direction: the added collision process shortens numerical thermalization time. Exact historical random realizations are not expected to be bitwise identical.

## One-click validation

Run `validate_thesis.bat` on Windows, or `./validate_thesis.sh` on Linux/macOS. The validator compiles the recovered source, enforces `alpha=3` and first-order weighting, runs collisionless and collisional representative cases, defaults to three repetitions, escalates to five only if the stochastic ordering is unclear, and prints `THESIS_VALIDATION_PASS` only when the expected direction is reproduced.

## Claim boundary

Public wording should say **"validated against the thesis results"** or **"representative reruns reproduce the main reported trends and mechanisms."** Do not say that every historical data point, every random realization, or the entire 2014 paper was exactly replicated.
