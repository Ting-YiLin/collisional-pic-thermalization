# Paper-to-code map

| Research element | Repository location |
|---|---|
| 1D electrostatic PIC particle/field update | `src/pic_2010_recovered.cpp` (`eq_of_motion`, density/field weighting functions) |
| Periodic particle boundary | `take_out_particle` |
| Large-angle / Turner collision | `MK_collision` |
| Pitch-angle electron-ion scattering | `PAS_collision` |
| Fermi-like / ball initial distribution | `ball_distribution` |
| Hermite-mode analysis | `Hermite_polynomials` |
| Multi-run averaging and linear fitting | `Ave_St`, `Lin_REG`, `main` |
| Historical parameter example | `config/paper_example_input.txt` |
| Modern end-to-end execution check | `scripts/smoke_test.py` |
| Published fitted scaling relations | `scripts/plot_reported_scalings.py` |

The 2014 paper later used a chi-square goodness-of-fit method for a more general determination of longitudinal and transverse thermalization times. That later analysis is documented in the paper but is not silently inserted into the recovered 2010 source.
