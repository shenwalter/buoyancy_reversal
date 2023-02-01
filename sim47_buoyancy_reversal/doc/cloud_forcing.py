#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

x1 = np.linspace(0, 1, 99)
z1 = np.linspace(0, 1, 99)

x, z = np.meshgrid(x1, z1)

psi = z * (z - 1)

niter = 100

# Forcing from Eq (13) of
# Grabowski 1993 "Cumulus entrainment, fine-scale mixing, and buoyancy reversal"
# https://doi.org/10.1002/qj.49711951305
pi = np.pi
hist_m = []
hist_diff = []
for m in range(niter):
    dpsi = (2 / pi)**3 * np.sin(
        (2 * m + 1) * pi * z) / (2 * m + 1)**3 * np.exp(-(2 * m + 1) * pi * x)
    diff = np.max(np.abs(dpsi))
    psi += dpsi
    hist_m.append(m)
    hist_diff.append(diff)

hist_m = np.array(hist_m)
hist_diff = np.array(hist_diff)
hist_errorbound = (2 / pi)**3 * 0.25 / (2 * hist_m + 1)**2
print("Last error bound: ", hist_errorbound[-1])

fig, ax = plt.subplots()
im = ax.imshow(np.flipud(psi), extent=(x1.min(), x1.max(), z1.min(), z1.max()))
ax.set_xticks([0, 1])
ax.set_yticks([0, 1])
ax.set_xlabel('x', labelpad=-8)
ax.set_ylabel('z', labelpad=-8)
cbar = fig.colorbar(im)
fig.savefig('psi.pdf')
plt.close(fig)

fig, ax = plt.subplots()
ax.plot(hist_m, hist_diff, label='diff')
ax.plot(hist_m, np.sum(hist_diff) - np.cumsum(hist_diff), label='error')
ax.plot(hist_m, hist_errorbound, label='errorbound')
ax.legend()
ax.set_xlabel('iter')
ax.set_yscale('log')
ax.set_ylim(1e-9, 1)
fig.savefig('hist.pdf')
plt.close(fig)
