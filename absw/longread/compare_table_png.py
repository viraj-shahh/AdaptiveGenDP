#!/usr/bin/env python3
"""Render a `golden compare` output file as a PNG table (for slides / chat screenshots).

usage: compare_table_png.py compare_<...>.txt [out.png]
Default output: same path with .png instead of .txt.
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

src = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else src.rsplit(".", 1)[0] + ".png"
lines = open(src).read().split("\n")
header = lines[0]                                    # "5000 pairs, mean qlen 1000, mean full-DP cells 1058226"
rows = {}                                            # (policy, W) -> (cells, % optimal)
for l in lines[2:]:
    p = l.split()
    if len(p) == 6:
        rows[(p[0], int(p[1]))] = (int(p[2]), p[3])
Ws = sorted({w for _, w in rows})
n_pairs = header.split()[0]
qlen = header.split("mean qlen ")[1].split(",")[0]
full = int(header.split("cells ")[1])

pols = ["rowband", "static", "adaptive"]
names = ["GenDP band\n(rowband)", "static", "adaptive"]
left = [[str(w)] + [rows[(p, w)][1] for p in pols] for w in Ws]
right = [[f"{rows[(p, w)][0] / 1000:.1f}k" for p in pols] for w in Ws]

fig, (axl, axr) = plt.subplots(1, 2, figsize=(9, 0.42 * (len(Ws) + 2.3)), dpi=200,
                               gridspec_kw={"width_ratios": [4.4, 3.6], "wspace": 0.03})


def draw(ax, data, labels, widths, title, adaptive_col):
    ax.axis("off")
    th = (len(data) + 1) / (len(data) + 1.7)      # table height; the rest is the merged group header
    t = ax.table(cellText=data, colLabels=labels, colWidths=widths, cellLoc="center", bbox=[0, 0, 1, th])
    t.auto_set_font_size(False); t.set_fontsize(10)
    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor("#d0d0d0")
        if r == 0:
            cell.set_facecolor("#2b4c7e"); cell.set_text_props(color="white", weight="bold")
        elif c == adaptive_col:
            cell.set_facecolor("#e3eefc"); cell.set_text_props(weight="bold")   # highlight adaptive
    ax.add_patch(plt.Rectangle((0, th), 1, 1 - th, transform=ax.transAxes, clip_on=False,
                               facecolor="#e6e6e6", edgecolor="#d0d0d0"))
    ax.text(0.5, (1 + th) / 2, title, transform=ax.transAxes, ha="center", va="center",
            fontsize=11, weight="bold")


draw(axl, left, ["W"] + names, [0.16, 0.32, 0.26, 0.26], "% optimal", 3)
draw(axr, right, names, [0.36, 0.32, 0.32], "cells/pair", 2)
fig.suptitle(f"{n_pairs} random $\\bf{{{qlen}}}$-base pieces of S. aureus ONT reads", fontsize=11, y=1.02)
fig.text(0.5, -0.02, f"full table = {full / 1e6:.2f}M cells/pair", ha="center", va="top",
         fontsize=9, color="#555")
fig.savefig(out, bbox_inches="tight", facecolor="white")
print(f"wrote {out}")
