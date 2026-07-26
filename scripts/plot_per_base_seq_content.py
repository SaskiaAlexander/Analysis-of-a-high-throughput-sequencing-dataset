#!/usr/bin/env python3
import glob
import os
import re

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FuncFormatter, MultipleLocator

# Folder containing the extracted *_fastqc folders produced by FastQC
FASTQC_DIR = r"C:\Users\alexa\OneDrive\Howest\BIT 01 Linux operating systems\SF\BIT09\data\02_fastqc"
# Folder for the output plot
OUT_DIR = r"C:\Users\alexa\OneDrive\Howest\BIT 01 Linux operating systems\SF\BIT09\plots"
# Nucleotide -> line color, matching FastQC's own convention (A=green, T=red, G=black, C=blue)
BASE_COLORS = [("A", "green"), ("T", "red"), ("G", "black"), ("C", "blue")]
BASES = [base for base, _ in BASE_COLORS]
Y_AXIS_PADDING = 5  # percentage points of headroom above/below the actual data range


def bin_to_pos(base_label: str) -> float:
    # FastQC bins later positions into ranges like "16-17" - use the midpoint so the
    # x-axis is a real bp scale instead of one categorical tick per row
    if "-" in base_label:
        start, end = base_label.split("-")
        return (float(start) + float(end)) / 2
    return float(base_label)


def parse_per_base_content(fastqc_dir: str) -> pd.DataFrame:
    # Each extracted FastQC folder contains a fastqc_data.txt with the raw per-base numbers
    with open(os.path.join(fastqc_dir, "fastqc_data.txt")) as f:
        text = f.read()

    # Extract only the "Per base sequence content" block from the report
    section = re.search(r">>Per base sequence content.*?\n(.*?)\n>>END_MODULE", text, re.S).group(1)

    # Skip the header line (starts with #) and split each row into columns
    rows = [line.split("\t") for line in section.strip().splitlines() if not line.startswith("#")]
    df = pd.DataFrame(rows, columns=["Base", "G", "A", "T", "C"])

    # Nucleotide percentages come in as text, convert to numeric for plotting
    df[["G", "A", "T", "C"]] = df[["G", "A", "T", "C"]].astype(float)
    df["Pos"] = df["Base"].apply(bin_to_pos)
    return df


def plot_sample(ax, df: pd.DataFrame, sample_name: str) -> None:
    # zorder=3 keeps the data lines drawn on top of the grid
    for base, color in BASE_COLORS:
        ax.plot(df["Pos"], df[base], label=base, color=color, zorder=3)

    # Numeric bp axis with a tick every 10 bp, matching MultiQC's "10 bp / 20 bp / ..." style
    ax.set_xlim(0, df["Pos"].max())
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x)} bp"))
    ax.tick_params(axis="x", labelsize=7)

    # Light grey gridlines like the MultiQC report, drawn behind the data lines
    ax.grid(True, color="lightgrey", linewidth=0.7, zorder=0)

    # Sample name identifies each row now that the axis titles are shared
    ax.set_title(sample_name, fontsize=9, loc="left")


# Collect all extracted FastQC result folders and parse each one's data once
fastqc_dirs = sorted(
    d for d in glob.glob(os.path.join(FASTQC_DIR, "*_fastqc")) if os.path.isdir(d)
)
samples = [
    (os.path.basename(d).replace("_fastqc", ""), parse_per_base_content(d))
    for d in fastqc_dirs
]

# Fit the y-axis to the actual data range instead of the full 0-100% scale
all_values = pd.concat(df[BASES] for _, df in samples)
y_min = max(0, all_values.values.min() - Y_AXIS_PADDING)
y_max = min(100, all_values.values.max() + Y_AXIS_PADDING)

# Arrange the samples in a grid of two columns and as many rows as needed
ncols = 2
nrows = (len(samples) + ncols - 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 2 * nrows), sharex=True, sharey=True, squeeze=False)
axes = axes.flatten()

for ax, (sample_name, df) in zip(axes, samples):
    plot_sample(ax, df, sample_name)
    ax.label_outer()  # only the left column keeps y labels and the bottom row keeps x labels

axes[0].set_ylim(y_min, y_max)  # sharey=True applies this to every subplot

# Reserve space around the subplot grid: left/bottom for the shared axis labels
fig.tight_layout(rect=[0.07, 0.09, 0.86, 0.86])

# Anchor the title, x-axis label and legend to the actual top/bottom edges of the
# subplot grid, instead of guessed constants, so they stay tight against the plots
top_edge = axes[0].get_position().y1
bottom_edge = axes[-1].get_position().y0

fig.suptitle("Per Base Sequence Content", y=top_edge + 0.09, fontsize=14)
fig.supxlabel("Position in read (bp)", y=bottom_edge - 0.10)
fig.supylabel("Sequence content (%)", x=0.055)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.855, top_edge))

# Make sure the output directory exists before saving into it
os.makedirs(OUT_DIR, exist_ok=True)
fig.savefig(os.path.join(OUT_DIR, "per_base_sequence_content_grid.png"), dpi=200, bbox_inches="tight")
