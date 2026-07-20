import pyphysics as phys
import matplotlib.pyplot as plt

plt.rcParams["axes.labelsize"] = 16

beams = ["7Li", "11Li"]
labels = [r"$^7Li$", r"$^{11}$Li"]
ebeams = [7 * 7.5, 11 * 7.5]

fig, ax = plt.subplots(1, 1, figsize=(4, 3.25), constrained_layout=True)
for i, beam in enumerate(beams):
    k = phys.Kinematics(f"{beam}(d,p)@{ebeams[i]}").get_line3()
    ax.plot(k[0], k[1], lw=1.25, label=f"{labels[i]}(d,p)")

# Axis settings
ax.set_xlabel(r"$\theta_{lab}$ [$\circ$]")
ax.set_ylabel(r"$E_{lab}$ [MeV]")
ax.set_ylim(0)
ax.set_xlim(0, 180)
ax.legend()

fig.savefig("./Outputs/kin_li.png", dpi=300)

plt.show()
