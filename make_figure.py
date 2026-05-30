"""
Publication-ready permeability figure for quercetin and its circulating
metabolites: TPSA vs cLogP, the two physicochemical axes that drive passive
membrane permeation (and underlie the SwissADME BOILED-Egg model).

Data are the real RDKit/SwissADME values computed for the manuscript. The shaded
band marks the passive-permeation-favourable region (TPSA <= 140 A^2, cLogP >= 0);
points are coloured by metabolite class so the aglycone-vs-conjugate split is
visually immediate.

Outputs a 300-dpi PNG and a vector PDF suitable for a supplementary figure.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Real values (RDKit; TPSA/cLogP as used in the manuscript table).
# name: (TPSA, cLogP, class, GI_absorption)
DATA = {
    "Quercetin\n(aglycone)":        (131.36, 1.99, "permeable", "High"),
    "Isorhamnetin\n(3'-O-methyl)":  (120.36, 2.29, "permeable", "High"),
    "Quercetin-3-\nO-sulfate":      (174.73, 1.46, "conjugate", "Low"),
    "Isoquercetin\n(3-O-glucoside)":(210.51, -0.54, "conjugate", "Low"),
    "Quercetin-3-\nO-glucuronide":  (227.58, -0.45, "conjugate", "Low"),
}

COLOURS = {"permeable": "#2a7f3f", "conjugate": "#b03a2e"}
LABELS = {"permeable": "Aglycone / methylated (passive permeation feasible)",
          "conjugate": "Conjugate (transporter-dependent)"}


def main():
    fig, ax = plt.subplots(figsize=(7.2, 5.4))

    # Shade the passive-permeation-favourable region: TPSA <= 140, cLogP >= 0
    xmax = 250
    ax.add_patch(Rectangle((0, 0), 140, 3.2, facecolor="#d8efdc",
                           edgecolor="none", zorder=0))
    ax.axvline(140, color="#2a7f3f", ls="--", lw=1, zorder=1)
    ax.text(141, -0.85, "TPSA = 140 \u00c5$^2$\n(Veber threshold)",
            fontsize=8, color="#2a7f3f", va="bottom")

    seen = set()
    for name, (tpsa, logp, cls, gi) in DATA.items():
        lbl = LABELS[cls] if cls not in seen else None
        seen.add(cls)
        ax.scatter(tpsa, logp, s=170, c=COLOURS[cls], edgecolor="black",
                   linewidth=0.8, zorder=3, label=lbl)
        # annotate with name + GI absorption
        dy = 0.16 if cls == "permeable" else -0.34
        ax.annotate(f"{name.strip()}\n(GI: {gi})", (tpsa, logp),
                    textcoords="offset points", xytext=(0, 10 if dy > 0 else -28),
                    ha="center", fontsize=7.2, zorder=4)

    ax.set_xlabel("Topological polar surface area, TPSA (\u00c5$^2$)", fontsize=11)
    ax.set_ylabel("Calculated lipophilicity (cLogP)", fontsize=11)
    ax.set_title("Passive permeation potential: quercetin aglycone vs circulating metabolites",
                 fontsize=11.5, pad=12)
    ax.axhline(0, color="grey", lw=0.7, ls=":", zorder=1)
    ax.set_xlim(90, xmax)
    ax.set_ylim(-1.2, 3.0)
    ax.legend(loc="upper right", fontsize=8.5, frameon=True)
    ax.grid(True, alpha=0.25, zorder=0)
    fig.tight_layout()

    fig.savefig("results/quercetin_permeability_figure.png", dpi=300)
    fig.savefig("results/quercetin_permeability_figure.pdf")
    print("Wrote results/quercetin_permeability_figure.png (300 dpi)")
    print("Wrote results/quercetin_permeability_figure.pdf (vector)")


if __name__ == "__main__":
    main()
