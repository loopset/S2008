import pyphysics as phys
import hist
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import uproot

file = uproot.open("../../Fits/20Mg/Outputs/preliminary_xs.root")
if file is None:
    raise ValueError("File is not reachable")

h = file["fHist"].to_hist()  # type: ignore

## Read good AZURE2 fit
azure = phys.parse_txt("../../Fits/20Mg/Azure/Outputs/AZUREOut_aa=1_R=1.out", ncols=9)
# And data
data = phys.parse_txt("../../Fits/20Mg/Azure/Inputs/lab_1425.dat", ncols=4)
# Transform to CM in direct kinematics
data[:, 0] *= 20.0 / 21
# And normalize according to azure output
data[:, 2:4] *= 7.8653740e-03


fig, ax = plt.subplots(figsize=(6, 5))
overflow = 250
# phys.utils.set_hist_overflow(h, overflow)
ret = h.plot(ax=ax, cmap="managua_r", cmax=overflow, cmin=1, rasterized=True)
ax.set_xlabel(r"$\theta_{lab}$ [$\circ$]")
ax.set_ylabel(r"$E_{beam,p}$ [MeV]")
ret[1].set_label(r"d$\sigma$/d$\Omega$ [mb/sr]")
ax.set_xlim(80)
ax.axvspan(xmin=140, xmax=145, color="crimson", alpha=0.25)

fig.tight_layout()
fig.savefig("./Outputs/xs.png", dpi=300)

# plt.close("all")
fig, ax = plt.subplots(figsize=(6, 4))
ax.errorbar(
    data[:, 0],
    data[:, 2],
    yerr=data[:, 3],
    ls="none",
    marker="s",
    ms=4,
    mfc="none",
    color="dodgerblue",
)
ax.plot(azure[:, 0], azure[:, 3], color="crimson", label="Azure2 fit")
ax.set_xlabel(r"$E_{CM}$ [MeV]")
ax.set_ylabel(r"$d\sigma/d\Omega$ [mb/sr]")

# Annotate
ax.annotate(
    r"$5/2^+$" + "\n" + r"$\Gamma = 48 eV$",
    xy=(1.2, 0.3),
    ha="center",
    va="center",
    fontsize=12,
)
ax.annotate(
    r"$1/2^+$" + "\n" + r"$\Gamma = 170 keV$",
    xy=(1.8, 0.3),
    ha="center",
    va="center",
    fontsize=12,
)

ax.legend()

fig.tight_layout()
fig.savefig("./Outputs/azure.png", dpi=300)

plt.show()
