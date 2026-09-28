"""Generate original visualizations of scaling relations reported in Lai et al. (2014).
This script does not digitize or copy journal figures. It plots the published fitted equations.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[1] / "assets"
OUT.mkdir(exist_ok=True)

# Fig. 3 text: omega_pe * tau_R = 2.42 * N_D^2 for collisionless 1D ES PIC.
nd = np.linspace(20, 200, 300)
tau = 2.42 * nd**2
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.plot(nd, tau, label=r"$\omega_{pe}\tau_R = 2.42N_D^2$")
ax.set_xlabel(r"Particles per Debye length, $N_D$")
ax.set_ylabel(r"Reported $\omega_{pe}\tau_R$")
ax.set_title("Collisionless scaling reported in Lai et al. (2014)")
ax.grid(alpha=0.25)
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "reported_collisionless_scaling.png", dpi=180)
plt.close(fig)

# Eqs. (14) and (16): ratio of transverse to longitudinal thermalization times.
x = np.logspace(-4, 4, 500)
large_angle = 1 + 0.2 * x**(-1.0)
pitch_angle = 1 + 0.33 * x**(-0.69)
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.loglog(x, large_angle, label="Large-angle collisions: Eq. (14)")
ax.loglog(x, pitch_angle, label="Pitch-angle scattering: Eq. (16)")
ax.axvline(1.0, linestyle="--", linewidth=1, label=r"$N_D^2\nu_0/\omega_{pe}=1$")
ax.set_xlabel(r"Scale parameter $N_D^2\nu_0/\omega_{pe}$")
ax.set_ylabel(r"$\tau_{R\perp}/\tau_{R\parallel}$")
ax.set_title("Published fitted scaling relations")
ax.grid(alpha=0.25, which="both")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "reported_collision_ratio_scalings.png", dpi=180)
plt.close(fig)
print("generated reported scaling figures")
