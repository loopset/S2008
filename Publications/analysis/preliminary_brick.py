import pyphysics as phys
import hist
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import uproot

file = uproot.open("../../Fits/20Mg/Outputs/preliminary_xs_cm.root")
if file is None:
    raise ValueError("File is not reachable")

h = file["fHist"].to_hist()  # type: ignore

# Read BRICK
brick = np.load("../../Fits/20Mg/Azure/Outputs/bayesian_results.npz")
print("Brick results: ")
for i, key in enumerate(brick["par_labels"]):
    print(
        f"{key} idx: {i}, val: {brick['par_medians'][i]} + {brick['par_ulows'][i]} - {brick['par_uups'][i]}"
    )

fig, ax = plt.subplots(figsize=(6, 4))
overflow = 500
# phys.utils.set_hist_overflow(h, overflow)
ret = h.plot(ax=ax, cmap="managua_r", cmax=overflow, cmin=1, rasterized=True)
ax.set_xlabel(r"$\theta_{CM}$ [$\circ$]")
ax.set_ylabel(r"$E_{CM}$ [MeV]")
ret[1].set_label(r"d$\sigma$/d$\Omega$ [mb/sr]")
ax.set_xlim(0)
# ax.axvspan(xmin=140, xmax=145, color="crimson", alpha=0.25)

fig.tight_layout()
fig.savefig("./Outputs/preliminary_brick_xs.png", dpi=300)


################################### AZURE plot
# plt.close("all")
fig, ax = plt.subplots(figsize=(6, 4))
ax.errorbar(
    brick["exp_ecm"],
    brick["exp_y"],
    yerr=brick["exp_uy"],
    ls="none",
    marker="s",
    ms=4,
    mfc="none",
    color="dodgerblue",
)
ax.plot(brick["exp_ecm"], brick["azure_y"], color="crimson", label="AZURE2")
ax.set_xlabel(r"$E_{CM}$ [MeV]")
ax.set_ylabel(r"$d\sigma/d\Omega$ [mb/sr]")

# Set angle
ax.set_title(r"$\theta_{CM}$ = 144.2$^{\circ}$")
# Separation energy
sep = -1.1
# Annotate
ax.annotate(
    rf"E(5/2$^+$) = {brick['par_medians'][3] - sep:.2f} MeV"
    + "\n"
    + rf"$\Gamma$ = {brick['par_medians'][4] * 1e-3:.2f} keV",
    xy=(1.1, 0.7),
    ha="center",
    va="center",
    fontsize=12,
)
ax.annotate(
    rf"E(1/2$^+$) = {brick['par_medians'][0] - sep:.2f} MeV"
    + "\n"
    + rf"$\Gamma$ = {brick['par_medians'][1] * 1e-3:.2f} keV",
    xy=(2, 0.5),
    ha="center",
    va="center",
    fontsize=12,
)

ax.legend()

fig.tight_layout()
fig.savefig("./Outputs/preliminary_brick_rmatrix.png", dpi=300)

plt.show()
