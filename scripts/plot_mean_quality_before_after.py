#!/usr/bin/env python3
################################################################################
# Author	: Saskia Alexander - Howest student
# Usage		: python plot_mean_quality_before_after.py
#
# Recreates the MultiQC "Mean Quality Scores" plot, overlaying the mean
# per-base quality before trimming (green) and after trimming (dark green) as
# two averaged lines in one figure. The per-position mean Phred values are read
# from the MultiQC export files (one value per sample per position) and averaged
# across the six samples for each condition.
################################################################################
# IMPORTS
################################################################################
import os
import re

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter
################################################################################
# CONFIGURATION
################################################################################
# Input data folder and output folder for the plot
DATA = r"/home/guest/Shared/BIT09/data"
OUT_DIR = r"/home/guest/Shared/BIT09/plots"

# (label, line colour, MultiQC export of the per-base sequence quality plot)
SERIES = [
    ("Before trimming", "#5cb85c", os.path.join(DATA, "02_fastqc", "fastqc_multiqc", "multiqc_data", "fastqc_per_base_sequence_quality_plot.txt")),
    ("After trimming",  "#1b5e20", os.path.join(DATA, "04_fastqc_trimmed", "fastqc_trimmed_multiqc", "multiqc_data", "fastqc_per_base_sequence_quality_plot.txt")),
]

# FastQC quality zones (low, high, colour) drawn as the shaded background
ZONES = [(28, 41, "#4caf50"), (20, 28, "#ff9800"), (15, 20, "#f44336")]
################################################################################
# FUNCTIONS
################################################################################
def mean_quality(path):
    # Each cell looks like "(12, 39.6766...)"; average the mean quality per
    # position across all samples into a single line.
    per_pos = {}
    with open(path) as f:
        next(f)  # skip the header row with the column indices
        for line in f:
            for cell in line.rstrip("\n").split("\t")[1:]:  # drop the sample name
                m = re.match(r"\(\s*(\d+)\s*,\s*([\d.]+)\s*\)", cell.strip())
                if m:
                    per_pos.setdefault(int(m.group(1)), []).append(float(m.group(2)))
    xs = sorted(per_pos)
    return xs, [sum(per_pos[x]) / len(per_pos[x]) for x in xs]
################################################################################
# BUILD THE PLOT
################################################################################
fig, ax = plt.subplots(figsize=(8.2, 4.2))

# Draw the shaded quality zones behind everything else
for low, high, colour in ZONES:
    ax.axhspan(low, high, color=colour, alpha=0.12, linewidth=0, zorder=0)

# Draw one averaged mean-quality line per condition
for label, colour, path in SERIES:
    x, y = mean_quality(path)
    ax.plot(x, y, color=colour, linewidth=2, label=label, zorder=3)
################################################################################
# AXES AND STYLING
################################################################################
ax.set_xlim(0, max(x))  # both series share the same positions
ax.set_ylim(15, 41)
ax.xaxis.set_major_locator(FixedLocator([20, 40, 60, 80]))  # no 0 bp / 100 bp label
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v)} bp"))
ax.set_axisbelow(True)
ax.grid(True, color="#d9d9d9", linewidth=0.8)
ax.set_xlabel("Position (bp)", color="#555555")
ax.set_ylabel("Phred Score", color="#555555")
ax.set_title("FastQC: Mean Quality Scores", color="#555555", fontsize=14, pad=18)
ax.tick_params(colors="#555555", length=0)  # grey tick labels, no tick marks
for spine in ax.spines.values():
    spine.set_visible(False)  # no frame around the plot
ax.legend(loc="lower right", labelcolor="#555555", frameon=True, edgecolor="#cccccc", facecolor="white", framealpha=0.9)
################################################################################
# SAVE OUTPUT
################################################################################
fig.tight_layout()
os.makedirs(OUT_DIR, exist_ok=True)
fig.savefig(os.path.join(OUT_DIR, "mean_quality_before_after.png"), dpi=200, bbox_inches="tight")
