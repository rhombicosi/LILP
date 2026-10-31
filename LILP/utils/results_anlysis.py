import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from constants_paths import *
import numpy as np
from scipy.optimize import curve_fit


# ============================================================
# Prepare results
# ============================================================
model= 'cb-xyz'
files_init = [
        f'{results_folder_name}/LILP_0_50_{model}.txt',
        f'{results_folder_name}/LILP_50_60_{model}.txt',
        f'{results_folder_name}/LILP_60_70_{model}.txt',
        f'{results_folder_name}/LILP_70_80_{model}.txt',
        f'{results_folder_name}/LILP_80_90_{model}.txt'
]

files_cleaned = [
    f'{results_folder_name}/LILP_0_50_{model}-cleaned.txt',
    f'{results_folder_name}/LILP_50_60_{model}-cleaned.txt',
    f'{results_folder_name}/LILP_60_70_{model}-cleaned.txt',
    f'{results_folder_name}/LILP_70_80_{model}-cleaned.txt',
    f'{results_folder_name}/LILP_80_90_{model}-cleaned.txt'
]

results_total = f'{results_folder_name}/LILP_0_90_{model}-total.txt'

def results_cleanup(results_file, results_clean):

    # Read the TXT file
    df = pd.read_csv(results_file, sep=r"\s+", engine="python")

    # 1st column width
    df["RNAseqname"] = df["RNAseqname"].str.replace(
        r"^(tRNA_tdbR\d+).*",
        r"\1",
        regex=True
    )

    df.to_csv(results_clean, sep="\t", index=False)
    df = pd.read_csv(results_clean, sep=r"\s+", engine="python")

    widths = {
        col: max(len(str(col)), df[col].astype(str).str.len().max())
        for col in df.columns
    }

    with open(results_clean, "w") as f:

        # Header
        f.write(
            " ".join(
                str(col).ljust(widths[col])
                for col in df.columns
            ).rstrip() + "\n"
        )

        # Rows
        for _, row in df.iterrows():
            values = []

            for col in df.columns:
                value = str(row[col])

                # First column = names, left aligned
                # Everything else = numbers, right aligned
                if col == df.columns[0]:
                    values.append(value.ljust(widths[col]))
                else:
                    values.append(value.rjust(widths[col]))

            f.write(" ".join(values).rstrip() + "\n")

    print(f"Saved: {results_clean}")

    # Check that it was read correctly
    print(df.head())
    print(df.columns)


for fi, fc in zip(files_init, files_cleaned):
    results_cleanup(fi, fc)

# ============================================================
# COMBINE FILES
# ============================================================
dfs = []

for file in files_cleaned:
    df = pd.read_csv(file, sep=r"\s+", engine="python")
    dfs.append(df)

# Combine all rows
combined = pd.concat(dfs, ignore_index=True)

# Save into one file
combined.to_csv(results_total, sep="\t", index=False)

print(f"Combined {len(combined)} sequences")
print(f"Saved: {results_total}")

df = pd.read_csv(results_total, sep=r"\s+", engine="python")

# ============================================================
# TIME + BEST
# ============================================================
    
def format_time(seconds):
        seconds = int(round(seconds))
        hours = seconds // 3600
        minutes = (seconds % 3600) // 60

        if hours > 0:
            if minutes > 0:
                return f"{hours}h{minutes}m"
            else:
                return f"{hours}h"

        return f"{minutes}m"

def TimeBest():

    # plt.rcParams["font.family"] = "Century Gothic"
    plt.rcParams["font.size"] = 24
    time_color = "#424862"
    best_color = "#FB9A63"
    
    # MEAN VALUES 
    mean_time = df.groupby("numofnts")["Time(s)"].mean()
    mean_best = df.groupby("numofnts")["Best(s)"].mean()

    # FIGURE
    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)

    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)

    # TIME — INDIVIDUAL VALUES
    ax.scatter(df["numofnts"], df["Time(s)"], color=time_color, alpha=0.30, s=30, label="Optimality proof")

    # TIME — MEAN
    ax.plot(mean_time.index, mean_time.values, color=time_color, linewidth=1.5, marker="o", markersize=2, label="Optimality proof mean")

    # =================================================
    # SHIFTED POWER-LAW FIT — TIME
    # =================================================

    x_time = mean_time.index.to_numpy(dtype=float)
    y_time = mean_time.values.astype(float)

    def shifted_power(x, a, n0, b, c):
        return a * np.maximum(x - n0, 1e-6) ** b + c

    # Initial parameters # a # n0 # b # c
    p0 = [ 1.0, x_time.min() - 5, 2.0, y_time.min()]

    # Bounds prevent pathological fits
    bounds = (
        [0, x_time.min() - 30, 0, 0],
        [np.inf, x_time.min() - 1e-6, 10, 4 * 3600]
    )

    params, _ = curve_fit(shifted_power, x_time, y_time, p0=p0, bounds=bounds, maxfev=100000)
    a, n0, b, c = params

    x_fit = np.linspace(x_time.min(), x_time.max(), 500)
    y_fit = shifted_power(x_fit, a, n0, b, c)

    ax.plot(x_fit, y_fit, color=time_color, linestyle="--", linewidth=1.5, label=f"Optimality proof trend")
    
    # =================================================
    # BEST INCUMBENT — INDIVIDUAL VALUES
    # =================================================
    ax.scatter(df["numofnts"], df["Best(s)"], color=best_color, alpha=0.30, s=30, label="Best incumbent")

    # =================================================
    # BEST INCUMBENT — MEAN
    # =================================================
    ax.plot(mean_best.index, mean_best.values, color=best_color, linewidth=1.5, marker="o", markersize=2, label="Best incumbent mean")

    # =================================================
    # SHIFTED POWER-LAW FIT — BEST
    # =================================================
    x_best = mean_best.index.to_numpy(dtype=float)
    y_best = mean_best.values.astype(float)

    p0 = [1.0, x_best.min() - 5, 2.0, y_best.min()]

    bounds = (
        [0, x_best.min() - 30, 0, 0],
        [np.inf, x_best.min() - 1e-6, 10, 4 * 3600]
    )

    params, _ = curve_fit(shifted_power, x_best, y_best, p0=p0, bounds=bounds, maxfev=100000)
    a, n0, b, c = params

    y_best_fit = shifted_power(x_fit, a, n0, b, c)

    ax.plot(x_fit, y_best_fit, color=best_color, linestyle="--", linewidth=1.5, label=f"Best incumbent trend")

    # =================================================
    # CUSTOM Y TICKS
    # =================================================
    interval = 30 * 60

    y_max = ax.get_ylim()[1]
    tick_max = np.ceil(y_max / interval) * interval

    ticks = np.arange(0, tick_max + interval, interval)

    ax.set_yticks([])

    for tick in ticks:

        if tick % 3600 == 0:
            tick_length = 0.025
            tick_width = 1.5
            fontsize = 14

            # Y-grid only at whole-hour positions
            ax.axhline(y=tick, alpha=0.25, linewidth=1)
        else:
            tick_length = 0.014
            tick_width = 1
            fontsize = 12

        ax.plot([-tick_length, 0], [tick, tick], transform=ax.get_yaxis_transform(), clip_on=False, color="black", linewidth=tick_width)

        ax.annotate(format_time(tick),
            xy=(-tick_length, tick),
            xycoords=ax.get_yaxis_transform(),
            xytext=(-10, 0),
            textcoords="offset points",
            rotation=45,
            rotation_mode="anchor",
            ha="right",
            va="bottom",
            fontsize=fontsize,
            annotation_clip=False
        )

    # =================================================
    # LAYOUT AXES LEGEND GRID SPINES
    # =================================================
    ax.set_xlabel("Sequence length (nucleotides)", fontsize=18, labelpad=15)
    ax.set_ylabel("Time", fontsize=18, labelpad=60)

    ax.set_title("Computational Scaling", fontsize=24, pad=20)

    ax.tick_params(axis="x", labelsize=12)

    ax.grid(True, alpha=0.25, linewidth=1)

    ax.legend(fontsize=14, frameon=False, loc="upper left")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.subplots_adjust(left=0.14, right=0.97, bottom=0.14, top=0.91)

    # =================================================
    # SAVE
    # =================================================
    plt.savefig("performance_vs_length.svg", dpi=100, bbox_inches="tight", transparent=True)

    plt.show()

# ============================================================
# MIP GAP — SEPARATE FIGURE
# ============================================================

def MipGap():

    gap_color = "#424862" 

    df["Gap_percent"] = df["Gap"] * 100

    mean_gap = ( df.groupby("numofnts")["Gap_percent"].mean().sort_index())

    # =================================================
    # FIGURE
    # =================================================

    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)

    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)

    # =================================================
    # MIP GAP — SCATTER
    # =================================================

    ax.scatter(df["numofnts"], df["Gap_percent"], color=gap_color, alpha=0.35, s = 30, label="MIP Gap")

    # =================================================
    # MIP GAP — MEAN
    # =================================================

    ax.plot( mean_gap.index, mean_gap.values, color=gap_color, linewidth=1.5, marker="o", markersize=2, label="Mean")

    # =================================================
    # TREND OF THE MEAN
    # =================================================

    x = mean_gap.index.to_numpy(dtype=float)
    y = mean_gap.values.astype(float)

    # Quadratic trend fitted to the mean values
    coefficients = np.polyfit(x, y, 2)
    polynomial = np.poly1d(coefficients)

    # Smooth x values for plotting
    x_smooth = np.linspace(x.min(), x.max(), 300)
    gap_trend = polynomial(x_smooth)

    # MIP Gap cannot be negative
    gap_trend = np.maximum(gap_trend, 0)

    # =================================================
    # TREND LINE 
    # =================================================

    ax.plot( x_smooth, gap_trend, color=gap_color, linestyle="--", linewidth=1.5, label="Trend")

    ax.set_xlabel( "Sequence length (nucleotides)", fontsize=18, labelpad=15)
    ax.set_ylabel("MIP Gap (%)", fontsize=18, labelpad=20)

    ax.set_title("Solution Quality Scaling", fontsize=24, pad=15)

    ax.tick_params(labelsize=12)

    ax.set_ylim(0, 50)

    ax.grid(True, alpha=0.25, linewidth=1)

    ax.legend(fontsize=14, frameon=False, loc="upper left")

    # =================================================
    # LAYOUT
    # =================================================
    plt.subplots_adjust(left=0.14, right=0.97, bottom=0.14, top=0.91)

    # =================================================
    # SAVE — TRANSPARENT
    # =================================================
    plt.savefig("mipgap_vs_length.svg", dpi=100, bbox_inches="tight", transparent=True)

    plt.show()

# ============================================================
# violin plot INF: NNILP, RNAstructure, RNAfold, UNAfold
# ============================================================
def violin_ILP_DP(results_total, imgfile):
    df_branch = pd.read_csv(f'{results_total}', sep=r"\s+")

    cols_branch = ["INFILP", "INFRNAstr", "INFRNAFold", "INFUNAFold"]
    data_branch = df_branch[cols_branch]
    data_branch = data_branch.rename(columns={'INFILP': 'NNILP', 'INFRNAstr': 'RNAstructure', 'INFRNAFold' : 'RNAFold', 'INFUNAFold': 'UNAFold'})
    means = data_branch.mean()

    # plt.rcParams["font.family"] = "Century Gothic"
    # plt.rcParams["font.size"] = 18
    colors = {
        'NNILP': '#BFC4D5',
        'RNAstructure': '#FEBC98',
        'RNAFold': '#424862',
        'UNAFold': '#FB9A63'
    }

    fig, ax = plt.subplots(figsize=(10, 7.5), dpi=100)
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)

    sns.violinplot(data=data_branch, palette=colors, inner='point')
    for i, m in enumerate(means):
        plt.text(i, m, f"{m:.3f}", color="black", ha="center", va="bottom", fontweight="bold", fontsize=14)

    plt.title("INF Distribution", fontsize=16)
    plt.ylabel("INF Value", fontsize=12)
    plt.xticks(rotation=-10, fontsize=12)
    plt.yticks(fontsize=12)

    plt.ylim(0, 1.0)
    plt.tight_layout()

    plt.savefig(imgfile, dpi=100, bbox_inches="tight", transparent=True)
    
    plt.show()

# calculate fails percentage
def failure_rate():
    methods = ["INFILP", "INFRNAstr", "INFRNAFold", "INFUNAFold"]

    # Sequences where ALL methods have INF = 0
    all_failed = (df[methods] == 0).all(axis=1)

    # Percentage of the complete dataset
    failure_percentage = all_failed.mean() * 100

    print(f"All-method failure rate: {failure_percentage:.1f}%")
    print(f"Number of sequences: {all_failed.sum()} / {len(df)}")


# ============================================================
#  execution
# ============================================================

# violin_ILP_DP(results_total, "inf_vs_length.svg")
# violin_ILP_DP(files_cleaned[0], "inf_vs_length_0_50.png")
# violin_ILP_DP(files_cleaned[1], "inf_vs_length_50_60.png")
# violin_ILP_DP(files_cleaned[2], "inf_vs_length_60_70.png")
# violin_ILP_DP(files_cleaned[3], "inf_vs_length_70_80.png")
# violin_ILP_DP(files_cleaned[4], "inf_vs_length_80_90.png")
violin_ILP_DP(f'{results_folder_name}/LILP_0_90_cb-xyz-total.txt', "inf_vs_length_0_90.svg")


TimeBest()
MipGap()

failure_rate()
