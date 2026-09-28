# LinkedIn copy

## Project title
Collisional PIC Numerical Thermalization — Restored & Validated Master's Research Code

## Project description
Recovered and validated the 1D electrostatic Particle-in-Cell simulation code I independently developed during my master's research in plasma physics. The repository preserves the chain from physical model to C++ implementation, numerical experiment, thesis results, and the later 2014 Physics of Plasmas publication. Representative reruns reproduce the main thesis trends and mechanisms; the thesis-validation preset uses alpha=3 as specified in the thesis.

## Post
I recently went back to a computational-physics project from my master's research and turned it into a public, reproducible archive.

The original project was a 1D electrostatic Particle-in-Cell simulation that I developed from scratch to study numerical thermalization and Monte-Carlo collision effects in Lorentz plasmas. That line of work later contributed to our 2014 Physics of Plasmas paper.

The restored C++ code now compiles in a modern environment and has been validated against the thesis results. Representative reruns reproduce the main reported trends and mechanisms, including the shortening of numerical thermalization time when pitch-angle collisions are added.

I kept the scope deliberately practical: this is not a new PIC framework or a claim that every historical random data point was exactly reproduced. It is a transparent record of the research chain — model → code → simulation → validation → publication.

Repository: {{REPO_URL}}
Paper DOI: https://doi.org/10.1063/1.4904307
Thesis: https://ir.lib.ncu.edu.tw/handle/987654321/25772
