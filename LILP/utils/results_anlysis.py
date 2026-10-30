import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from constants_paths import *
import numpy as np
from matplotlib.ticker import FuncFormatter
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.interpolate import make_interp_spline
from scipy.signal import savgol_filter


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

model= 'cb-xyz'
files_init = [
        f'{results_folder_name}/LILP_0_50_{model}.txt',
        f'{results_folder_name}/LILP_50_60_{model}.txt',
        # f'{results_folder_name}/LILP_60_70_{model}.txt',
        # f'{results_folder_name}/LILP_70_80_{model}.txt',
        # f'{results_folder_name}/LILP_80_90_{model}.txt'
]

files_cleaned = [
    f'{results_folder_name}/LILP_0_50_{model}-cleaned.txt',
    f'{results_folder_name}/LILP_50_60_{model}-cleaned.txt',
    # f'{results_folder_name}/LILP_60_70_{model}-cleaned.txt',
    # f'{results_folder_name}/LILP_70_80_{model}-cleaned.txt',
    # f'{results_folder_name}/LILP_80_90_{model}-cleaned.txt'
]

results_total = f'{results_folder_name}/LILP_0_60_{model}-total.txt'

for fi, fc in zip(files_init, files_cleaned):
    results_cleanup(fi, fc)

### COMBINE FILES

dfs = []

for file in files_cleaned:
    df = pd.read_csv(file, sep=r"\s+", engine="python")
    dfs.append(df)

# Combine all rows
combined = pd.concat(dfs, ignore_index=True)

# Save one long file
combined.to_csv(
    results_total,
    sep="\t",
    index=False
)

print(f"Combined {len(combined)} sequences")
print(f"Saved: {results_total}")

df = pd.read_csv(results_total, sep=r"\s+", engine="python")

# methods = ["INFILP", "INFRNAstr", "INFRNAFold", "INFUNAFold"]

# # Sequences where ALL methods have INF = 0
# all_failed = (df[methods] == 0).all(axis=1)

# # Percentage of the complete dataset
# failure_percentage = all_failed.mean() * 100

# print(f"All-method failure rate: {failure_percentage:.1f}%")
# print(f"Number of sequences: {all_failed.sum()} / {len(df)}")

# Create 3 scatter plot graphs: Time, Best Incumbent, MIPGap
def create_scatter_plots():
    fig, axes = plt.subplots(3, 1, figsize=(9, 12))

    # Time
    axes[0].scatter(df["numofnts"], df["Time(s)"])
    axes[0].set_xlabel("Number of nucleotides")
    axes[0].set_ylabel("Time (s)")
    axes[0].set_title("Time vs Number of Nucleotides")
    axes[0].grid(True)

    # Best
    axes[1].scatter(df["numofnts"], df["Best(s)"])
    axes[1].set_xlabel("Number of nucleotides")
    axes[1].set_ylabel("Best (s)")
    axes[1].set_title("Best Time vs Number of Nucleotides")
    axes[1].grid(True)

    # Gap
    axes[2].scatter(df["numofnts"], df["Gap"])
    axes[2].set_xlabel("Number of nucleotides")
    axes[2].set_ylabel("Gap")
    axes[2].set_title("Gap vs Number of Nucleotides")
    axes[2].grid(True)

    plt.tight_layout()
    plt.savefig("performance_vs_numofnts.png", dpi=300)
    plt.show()

### SCATTER PLOT + LINEAR TREND
def scatter_with_average_line(x, y, xlabel, ylabel, title, filename):
    plt.figure(figsize=(9, 6))

    # Scatter points
    plt.scatter(x, y, label="Measurements")

    # Linear fit
    coefficients = np.polyfit(x, y, 1)
    linear_fit = np.poly1d(coefficients)

    x_line = np.linspace(x.min(), x.max(), 200)
    y_line = linear_fit(x_line)

    # Trend line
    plt.plot(
        x_line,
        y_line,
        linewidth=2,
        label="Average trend"
    )

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()

# # Time
# scatter_with_average_line(
#     df["numofnts"],
#     df["Time(s)"],
#     "Number of nucleotides",
#     "Time (s)",
#     "Time vs Number of Nucleotides",
#     "time_vs_numofnts.png"
# )

# ============================================================
# FUNCTION FOR TIME / BEST
# ============================================================

def plot_time_metric(df, y_column, y_label, title, filename):
    
    x = df["numofnts"]
    y = df[y_column]

    plt.figure(figsize=(10, 7))

    # --------------------------------------------------------
    # 1. Individual measurements
    # --------------------------------------------------------
    plt.scatter(
        x,
        y,
        color="steelblue",
        alpha=0.45,
        s=45,
        label="Measurements"
    )

    # --------------------------------------------------------
    # 2. Mean for each number of nucleotides
    # --------------------------------------------------------
    mean_values = (
        df.groupby("numofnts")[y_column]
          .mean()
          .sort_index()
    )

    plt.plot(
        mean_values.index,
        mean_values.values,
        color="darkorange",
        marker="o",
        markersize=6,
        linewidth=2.5,
        label="Mean per nucleotide count"
    )

    # --------------------------------------------------------
    # 3. Nonlinear trend line
    # --------------------------------------------------------
    coefficients = np.polyfit(x, y, 2)
    polynomial = np.poly1d(coefficients)

    x_trend = np.linspace(x.min(), x.max(), 300)
    y_trend = polynomial(x_trend)

    plt.plot(
        x_trend,
        y_trend,
        color="crimson",
        linestyle="--",
        linewidth=2.5,
        label="Nonlinear trend"
    )

    # --------------------------------------------------------
    # 4. Formatting
    # --------------------------------------------------------
    plt.xlabel("Number of nucleotides", fontsize=12)
    plt.ylabel(y_label, fontsize=12)

    plt.title(
        title,
        fontsize=14,
        fontweight="bold"
    )

    plt.grid(
        True,
        alpha=0.25
    )

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        filename,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


# ============================================================
# FUNCTION FOR MIP GAP
# ============================================================

def plot_mipgap(df):
    
    x = df["numofnts"]
    y = df["Gap"] * 100     # convert 0–1 to 0–100 %

    plt.figure(figsize=(10, 7))

    # --------------------------------------------------------
    # 1. Individual measurements
    # --------------------------------------------------------
    plt.scatter(
        x,
        y,
        color="steelblue",
        alpha=0.45,
        s=45,
        label="Measurements"
    )

    # --------------------------------------------------------
    # 2. Mean MIP Gap for each number of nucleotides
    # --------------------------------------------------------
    mean_gap = (
        df.assign(Gap_percent=df["Gap"] * 100)
          .groupby("numofnts")["Gap_percent"]
          .mean()
          .sort_index()
    )

    plt.plot(
        mean_gap.index,
        mean_gap.values,
        color="darkorange",
        marker="o",
        markersize=6,
        linewidth=2.5,
        label="Mean per nucleotide count"
    )

    # --------------------------------------------------------
    # 3. Nonlinear trend line
    # --------------------------------------------------------
    coefficients = np.polyfit(x, y, 2)
    polynomial = np.poly1d(coefficients)

    x_trend = np.linspace(x.min(), x.max(), 300)
    y_trend = polynomial(x_trend)

    plt.plot(
        x_trend,
        y_trend,
        color="crimson",
        linestyle="--",
        linewidth=2.5,
        label="Nonlinear trend"
    )

    # --------------------------------------------------------
    # 4. Formatting
    # --------------------------------------------------------
    plt.xlabel(
        "Number of nucleotides",
        fontsize=12
    )

    plt.ylabel(
        "MIP Gap (%)",
        fontsize=12
    )

    plt.title(
        "MIP Gap vs Number of Nucleotides",
        fontsize=14,
        fontweight="bold"
    )

    plt.ylim(-2, 102)

    plt.yticks(
        np.arange(0, 101, 10)
    )

    plt.grid(
        True,
        alpha=0.25
    )

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        "mipgap_vs_numofnts.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


# ============================================================
# CREATE THE PLOTS
# ============================================================

# # TIME
# plot_time_metric(
#     df,
#     "Time(s)",
#     "Time (s)",
#     "Time vs Number of Nucleotides",
#     "time_vs_numofnts.png"
# )
# # BEST
# plot_time_metric(
#     df,
#     "Best(s)",
#     "Best Time (s)",
#     "Best Time vs Number of Nucleotides",
#     "best_vs_numofnts.png"
# )
# # MIP GAP
# plot_mipgap(df)


# ============================================================
# TIME + BEST
# ============================================================

    
def format_time(seconds, pos):
    """Convert seconds to h/m/s format without spaces."""
    seconds = int(round(seconds))

    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    secs = seconds % 60

    if hours > 0:
        return f"{hours}h{minutes}m{secs}s"
    elif minutes > 0:
        return f"{minutes}m{secs}s"
    else:
        return f"{secs}s"

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

# def TimeBest():

#     # =================================================
#     # COLORS
#     # =================================================

#     time_color = "tab:blue"
#     best_color = "tab:orange"

#     # =================================================
#     # MEAN VALUES
#     # =================================================

#     mean_time = df.groupby("numofnts")["Time(s)"].mean()
#     mean_best = df.groupby("numofnts")["Best(s)"].mean()

#     # =================================================
#     # FIGURE
#     # =================================================

#     fig, ax = plt.subplots(figsize=(12, 7))

#     fig.patch.set_alpha(0)
#     ax.patch.set_alpha(0)

#     # =================================================
#     # TIME
#     # =================================================

#     ax.scatter(
#         df["numofnts"],
#         df["Time(s)"],
#         color=time_color,
#         alpha=0.35,
#         label="Optimality proof"
#     )

#     ax.plot(
#         mean_time.index,
#         mean_time.values,
#         color=time_color,
#         linewidth=2.5,
#         label="Mean Optimality proof"
#     )

#     # Power-law regression:
#     # T(n) = a * n^b

#     log_time_x = np.log(mean_time.index)
#     log_time_y = np.log(mean_time.values)

#     time_result = linregress(
#         log_time_x,
#         log_time_y
#     )

#     time_b = time_result.slope
#     time_a = np.exp(time_result.intercept)
#     time_r2 = time_result.rvalue ** 2

#     x_smooth = np.linspace(
#         mean_time.index.min(),
#         mean_time.index.max(),
#         300
#     )

#     time_trend = time_a * x_smooth ** time_b

#     ax.plot(
#         x_smooth,
#         time_trend,
#         color=time_color,
#         linestyle="--",
#         linewidth=3,
#         label=f"Power-law Optimality proof ($b={time_b:.2f}$, $R^2={time_r2:.3f}$)"
#     )

#     # =================================================
#     # BEST
#     # =================================================

#     ax.scatter(
#         df["numofnts"],
#         df["Best(s)"],
#         color=best_color,
#         alpha=0.35,
#         label="Best incumbent"
#     )

#     ax.plot(
#         mean_best.index,
#         mean_best.values,
#         color=best_color,
#         linewidth=2.5,
#         label="Mean Best incumbent"
#     )

#     # Power-law regression:
#     # T(n) = a * n^b

#     log_best_x = np.log(mean_best.index)
#     log_best_y = np.log(mean_best.values)

#     best_result = linregress(
#         log_best_x,
#         log_best_y
#     )

#     best_b = best_result.slope
#     best_a = np.exp(best_result.intercept)
#     best_r2 = best_result.rvalue ** 2

#     best_trend = best_a * x_smooth ** best_b

#     ax.plot(
#         x_smooth,
#         best_trend,
#         color=best_color,
#         linestyle="--",
#         linewidth=3,
#         label=f"Power-law Best incumbent ($b={best_b:.2f}$, $R^2={best_r2:.3f}$)"
#     )

#     # =================================================
#     # Y-AXIS: TIME TICKS
#     # =================================================

#     # Get y-axis limits
#     y_min, y_max = ax.get_ylim()

#     # 30-minute intervals
#     interval = 30 * 60

#     first_tick = np.ceil(y_min / interval) * interval

#     ticks = np.arange(
#         first_tick,
#         y_max + interval,
#         interval
#     )

#     # Remove normal y-axis ticks
#     ax.set_yticks([])

#     # -------------------------------------------------
#     # Draw custom ticks
#     # -------------------------------------------------

#     for tick in ticks:

#         # 1 hour = long tick
#         if tick % 3600 == 0:
#             tick_length = 0.025
#             tick_width = 1.6

#         # 30 minutes = shorter tick
#         else:
#             tick_length = 0.014
#             tick_width = 1.1

#         ax.plot(
#             [-tick_length, 0],
#             [tick, tick],
#             transform=ax.get_yaxis_transform(),
#             clip_on=False,
#             color="black",
#             linewidth=tick_width
#         )

#         # -------------------------------------------------
#         # Label
#         #
#         # Bottom/right end of text starts at tick height.
#         # -------------------------------------------------

#         ax.annotate(
#             format_time(tick),
#             xy=(-tick_length, tick),
#             xycoords=ax.get_yaxis_transform(),
#             xytext=(-8, 0),
#             textcoords="offset points",
#             rotation=45,
#             rotation_mode="anchor",
#             ha="right",
#             va="bottom",
#             fontsize=10,
#             annotation_clip=False
#         )

#     # =================================================
#     # AXIS LABELS
#     # =================================================

#     ax.set_xlabel(
#         "Sequence length (nucleotides)"
#     )

#     ax.set_ylabel(
#         "Optimization time (s)",
#         labelpad=50
#     )

#     ax.set_title(
#         "Optimization Time vs. Sequence Length"
#     )

#     # =================================================
#     # GRID
#     # =================================================

#     ax.grid(
#         True,
#         alpha=0.3
#     )

#     # =================================================
#     # LEGEND
#     # =================================================

#     ax.legend()

#     # =================================================
#     # LAYOUT
#     # =================================================

#     # Extra space on the left so the tick labels
#     # cannot overlap the y-axis title.
#     plt.subplots_adjust(
#         left=0.25,
#         right=0.97,
#         bottom=0.12,
#         top=0.94
#     )

#     # =================================================
#     # SAVE
#     # =================================================

#     plt.savefig(
#         "performance_vs_length.png",
#         dpi=300,
#         bbox_inches="tight",
#         transparent=True
#     )

#     plt.show()

#     # =================================================
#     # STATISTICAL RESULTS
#     # =================================================

#     print("\nPower-law scaling results")
#     print("--------------------------------")

#     print(
#         f"Time: T(n) = {time_a:.4g} × n^{time_b:.4f}"
#     )

#     print(
#         f"Time exponent (b): {time_b:.4f}"
#     )

#     print(
#         f"Time R²: {time_r2:.4f}"
#     )

#     print()

#     print(
#         f"Best: T(n) = {best_a:.4g} × n^{best_b:.4f}"
#     )

#     print(
#         f"Best exponent (b): {best_b:.4f}"
#     )

#     print(
#         f"Best R²: {best_r2:.4f}"
#     )



def TimeBest():

    plt.rcParams["font.family"] = "Century Gothic"
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
    ax.scatter(
        df["numofnts"],
        df["Time(s)"],
        color=time_color,
        alpha=0.30,
        s=30,
        label="Optimality proof"
    )

    # TIME — MEAN
    ax.plot(
        mean_time.index,
        mean_time.values,
        color=time_color,
        linewidth=1.5,
        marker="o",
        markersize=2,
        label="Optimality proof mean"
    )

    # =================================================
    # TIME — POWER-LAW SCALING FIT
    # =================================================

    # x_time = mean_time.index.to_numpy()
    # y_time = mean_time.values

    # result_time = linregress(
    #     np.log(x_time),
    #     np.log(y_time)
    # )

    # time_b = result_time.slope
    # time_a = np.exp(result_time.intercept)
    # time_r2 = result_time.rvalue ** 2

    # x_fit = np.linspace(
    #     x_time.min(),
    #     x_time.max(),
    #     500
    # )

    # y_fit = time_a * x_fit ** time_b

    # # Show only the part within the experimental time limit
    # time_limit = 4 * 60 * 60

    # valid = y_fit <= time_limit

    # ax.plot(
    #     x_fit[valid],
    #     y_fit[valid],
    #     color=time_color,
    #     linestyle="--",
    #     linewidth=3.5,
    #     label=f"Scaling fit ($b={time_b:.2f}$, $R^2={time_r2:.3f}$)"
    # )



    # =================================================
    # SHIFTED POWER-LAW FIT — TIME
    # =================================================

    x_time = mean_time.index.to_numpy(dtype=float)
    y_time = mean_time.values.astype(float)

    def shifted_power(x, a, n0, b, c):
        return a * np.maximum(x - n0, 1e-6) ** b + c

    # Initial parameters
    p0 = [
        1.0,                  # a
        x_time.min() - 5,     # n0
        2.0,                  # b
        y_time.min()          # c
    ]

    # Bounds prevent pathological fits
    bounds = (
        [0, x_time.min() - 30, 0, 0],
        [np.inf, x_time.min() - 1e-6, 10, 4 * 3600]
    )

    params, _ = curve_fit(
        shifted_power,
        x_time,
        y_time,
        p0=p0,
        bounds=bounds,
        maxfev=100000
    )

    a, n0, b, c = params

    x_fit = np.linspace(
        x_time.min(),
        x_time.max(),
        500
    )

    y_fit = shifted_power(
        x_fit,
        a, n0, b, c
    )

    ax.plot(
        x_fit,
        y_fit,
        color=time_color,
        linestyle="--",
        linewidth=1.5,
        label=f"Optimality proof trend"
    )
    
    # =================================================
    # BEST INCUMBENT — INDIVIDUAL VALUES
    # =================================================

    ax.scatter(
        df["numofnts"],
        df["Best(s)"],
        color=best_color,
        alpha=0.30,
        s=30,
        label="Best incumbent"
    )

    # =================================================
    # BEST INCUMBENT — MEAN
    # =================================================

    ax.plot(
        mean_best.index,
        mean_best.values,
        color=best_color,
        linewidth=1.5,
        marker="o",
        markersize=2,
        label="Best incumbent mean"
    )

    # =================================================
    # BEST INCUMBENT — POWER-LAW SCALING FIT
    # =================================================

    # x_best = mean_best.index.to_numpy()
    # y_best = mean_best.values

    # result_best = linregress(
    #     np.log(x_best),
    #     np.log(y_best)
    # )

    # best_b = result_best.slope
    # best_a = np.exp(result_best.intercept)
    # best_r2 = result_best.rvalue ** 2

    # y_best_fit = best_a * x_fit ** best_b

    # valid_best = y_best_fit <= time_limit

    # ax.plot(
    #     x_fit[valid_best],
    #     y_best_fit[valid_best],
    #     color=best_color,
    #     linestyle="--",
    #     linewidth=3.5,
    #     label=f"Scaling fit ($b={best_b:.2f}$, $R^2={best_r2:.3f}$)"
    # )

    # =================================================
    # SHIFTED POWER-LAW FIT — BEST
    # =================================================

    x_best = mean_best.index.to_numpy(dtype=float)
    y_best = mean_best.values.astype(float)

    p0 = [
        1.0,
        x_best.min() - 5,
        2.0,
        y_best.min()
    ]

    bounds = (
        [0, x_best.min() - 30, 0, 0],
        [np.inf, x_best.min() - 1e-6, 10, 4 * 3600]
    )

    params, _ = curve_fit(
        shifted_power,
        x_best,
        y_best,
        p0=p0,
        bounds=bounds,
        maxfev=100000
    )

    a, n0, b, c = params

    y_best_fit = shifted_power(
        x_fit,
        a, n0, b, c
    )

    ax.plot(
        x_fit,
        y_best_fit,
        color=best_color,
        linestyle="--",
        linewidth=1.5,
        label=f"Best incumbent trend"
    )



    # =================================================
    # CUSTOM Y TICKS
    # =================================================

    interval = 30 * 60

    y_max = ax.get_ylim()[1]
    tick_max = np.ceil(y_max / interval) * interval

    ticks = np.arange(
        0,
        tick_max + interval,
        interval
    )

    ax.set_yticks([])

    for tick in ticks:

        if tick % 3600 == 0:
            tick_length = 0.025
            tick_width = 1.5
            fontsize = 14

            # Y-grid only at whole-hour positions
            ax.axhline(
                y=tick,
                alpha=0.25,
                linewidth=1
            )
        else:
            tick_length = 0.014
            tick_width = 1
            fontsize = 12

        ax.plot(
            [-tick_length, 0],
            [tick, tick],
            transform=ax.get_yaxis_transform(),
            clip_on=False,
            color="black",
            linewidth=tick_width
        )

        ax.annotate(
            format_time(tick),
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
    # AXES
    # =================================================

    ax.set_xlabel(
        "Sequence length (nucleotides)",
        fontsize=18,
        labelpad=15
    )

    ax.set_ylabel(
        "Time",
        fontsize=18,
        labelpad=60
    )

    ax.set_title(
        "Computational Scaling",
        fontsize=24,
        pad=20
    )

    ax.tick_params(
        axis="x",
        labelsize=12
    )

    # =================================================
    # GRID
    # =================================================

    ax.grid(
        True,
        alpha=0.25,
        linewidth=1
    )

    # =================================================
    # LEGEND — TRANSPARENT
    # =================================================

    ax.legend(
        fontsize=14,
        frameon=False,
        loc="upper left"
    )

    # =================================================
    # SPINES
    # =================================================

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # =================================================
    # LAYOUT
    # =================================================

    plt.subplots_adjust(
        left=0.14,
        right=0.97,
        bottom=0.14,
        top=0.91
    )

    # =================================================
    # SAVE
    # =================================================

    plt.savefig(
        "performance_vs_length.svg",
        dpi=100,
        bbox_inches="tight",
        transparent=True
    )

    plt.show()

# ============================================================
# MIP GAP — SEPARATE FIGURE
# ============================================================

def MipGap():

    gap_color = "#424862" 

    df["Gap_percent"] = df["Gap"] * 100

    mean_gap = (
        df.groupby("numofnts")["Gap_percent"]
        .mean()
        .sort_index()
    )

    # =================================================
    # FIGURE
    # =================================================

    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)

    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)

    # =================================================
    # MIP GAP — SCATTER
    # =================================================

    ax.scatter(
        df["numofnts"],
        df["Gap_percent"],
        color=gap_color,
        alpha=0.35,
        s = 30,
        label="MIP Gap"
    )

    # =================================================
    # MIP GAP — MEAN
    # =================================================

    ax.plot(
        mean_gap.index,
        mean_gap.values,
        color=gap_color,
        linewidth=1.5,
        marker="o",
        markersize=2,
        label="Mean"
    )

    # =================================================
    # TREND OF THE MEAN
    #
    # Uses ALL mean values, including Gap = 0.
    # No values are removed.
    # =================================================

    x = mean_gap.index.to_numpy(dtype=float)
    y = mean_gap.values.astype(float)

    # Quadratic trend fitted to the mean values
    coefficients = np.polyfit(
        x,
        y,
        2
    )

    polynomial = np.poly1d(
        coefficients
    )

    # Smooth x values for plotting
    x_smooth = np.linspace(
        x.min(),
        x.max(),
        300
    )

    gap_trend = polynomial(
        x_smooth
    )

    # MIP Gap cannot be negative
    gap_trend = np.maximum(
        gap_trend,
        0
    )

    # =================================================
    # TREND LINE
    # =================================================

    ax.plot(
        x_smooth,
        gap_trend,
        color=gap_color,
        linestyle="--",
        linewidth=1.5,
        label="Trend"
    )

    # =================================================
    # AXIS LABELS
    # =================================================

    ax.set_xlabel(
        "Sequence length (nucleotides)",
        fontsize=18,
        labelpad=15
    )

    ax.set_ylabel(
        "MIP Gap (%)",
        fontsize=18,
        labelpad=20
    )

    ax.set_title(
        "Solution Quality Scaling",
        fontsize=24,
        pad=15
    )

    ax.tick_params(
        labelsize=12
    )

    # =================================================
    # Y AXIS
    # =================================================

    ax.set_ylim(
        0,
        50
    )

    # =================================================
    # GRID
    # =================================================

    ax.grid(
        True,
        alpha=0.25,
        linewidth=1
    )

    # =================================================
    # LEGEND
    # =================================================

    ax.legend(
        fontsize=14,
        frameon=False,
        loc="upper left"
    )

    # =================================================
    # LAYOUT
    # =================================================

    plt.subplots_adjust(
        left=0.14,
        right=0.97,
        bottom=0.14,
        top=0.91
    )

    # =================================================
    # SAVE — TRANSPARENT
    # =================================================

    plt.savefig(
        "mipgap_vs_length.svg",
        dpi=100,
        bbox_inches="tight",
        transparent=True
    )

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

failure_rate()

# create violin plot INF of BLILP RNAstructure, RNAfold, UNAfold

def violin_ILP_DP(results_total, imgfile):
    df_branch = pd.read_csv(f'{results_total}', sep=r"\s+")

    cols_branch = ["INFILP", "INFRNAstr", "INFRNAFold", "INFUNAFold"]
    data_branch = df_branch[cols_branch]
    data_branch = data_branch.rename(columns={'INFILP': 'NNILP', 'INFRNAstr': 'RNAstructure', 'INFRNAFold' : 'RNAFold', 'INFUNAFold': 'UNAFold'})
    means = data_branch.mean()

    plt.rcParams["font.family"] = "Century Gothic"
    # plt.rcParams["font.size"] = 18
    colors = {
        'NNILP': '#BFC4D5',
        'RNAstructure': '#FEBC98',
        'RNAFold': '#424862',
        'UNAFold': '#FB9A63'
    }

    # fig, ax = plt.subplots(figsize=(10, 7.5), dpi=100)
    fig, ax = plt.subplots(figsize=(4.8, 3.6), dpi=100)
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

    plt.savefig(
        imgfile,
        dpi=100,
        bbox_inches="tight",
        transparent=True
    )
    
    plt.show()



    # df_melted = df.melt(id_vars=["RNAseqname"], value_vars=["F1ILP", "F1RNAstr", "F1UNAFold"],
    #                     var_name="Metric", value_name="Value")

    # sns.violinplot(x="Metric", y="Value", data=df_melted)
    # plt.title("RNA Structure Prediction Scores")
    # plt.ylim(0, 1.0)
    # plt.show()

# violin_ILP_DP(results_total, "inf_vs_length.svg")
violin_ILP_DP(files_cleaned[0], "inf_vs_length_0_50.png")
violin_ILP_DP(files_cleaned[1], "inf_vs_length_50_60.png")
# violin_ILP_DP(files_cleaned[2], "inf_vs_length_60_70.png")
# violin_ILP_DP(files_cleaned[3], "inf_vs_length_70_80.png")
violin_ILP_DP(f'{results_folder_name}/LILP_0_60_cb-xyz-total.txt', "inf_vs_length_0_60.svg")


TimeBest()
MipGap()
