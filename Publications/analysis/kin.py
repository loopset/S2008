import pyphysics as phys
import hist
import numpy as np
import matplotlib.pyplot as plt
import uproot
import matplotlib.axes as mplaxes

# Only front so far
data = uproot.open(
    "../../PostAnalysis/Outputs/tree_ex_20Mg_p_p_front.root:Final_Tree"
).arrays(["fThetaLight", "EVertex", "IsRANSAC"])  # type: ignore

hs = []
for i, ransac in enumerate([None, True]):
    if ransac is not None:
        mask = data["IsRANSAC"] == ransac
    else:
        mask = [True] * len(data)
    h = (
        hist.Hist.new.Reg(180, 0, 90, label=r"$\theta_{lab}$ [$\circ$]")
        .Reg(200, 0, 20, label=r"$E_{lab}$ [MeV]")
        .Double()
    )
    h.fill(data[mask]["fThetaLight"], data[mask]["EVertex"])
    hs.append(h)

fig, ax = plt.subplots(figsize=(6, 5))
ax: mplaxes.Axes
for i, h in enumerate(hs):
    cmap = "Reds_r" if i == 1 else "managua_r"
    alpha = 0.75 if i == 1 else 1
    h.plot(
        ax=ax,
        cmin=1,
        cbar=True if i == 0 else False,
        flow=False,
        cmap=cmap,
        alpha=alpha,
        rasterized=True,
    )


# plot theoretical line
theo = phys.Kinematics("20Mg(p,p)@84.85").get_line3()
ax.plot(theo[0], theo[1])

ax.set_title("Front events")

# annotations
ax.annotate(
    "Continuity",
    xy=(30, 2),
    xytext=(60, 1.5),
    ha="center",
    va="center",
    fontsize=14,
    color="dodgerblue",
    arrowprops=dict(arrowstyle="->"),
)

ax.annotate(
    "RANSAC",
    xy=(18, 11),
    xytext=(40, 13),
    ha="center",
    va="center",
    fontsize=14,
    color="crimson",
    arrowprops=dict(arrowstyle="->"),
)


fig.tight_layout()

fig.savefig("./Outputs/kin.png", dpi=300)
plt.show()
