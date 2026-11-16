#!/usr/bin/env python3
"""
README_CODE_BLOCK:
Purpose:
  Make a paper-style SIMP exclusion plot using the hps-simpy external-limit
  data files, with optional HPS reach/exclusion contours overlaid.

  The external/reference part of the plot follows the limits.py-style figure:
    - BaBar, Orsay, and E137 are borderless gray excluded regions.
    - CMB is a light gray solid contour line.
    - m_V/m_pi = 1.8 is a thick black line.
    - m_V/m_pi = 1.6 is a dashed gray line.
    - The top labels follow the HPS paper-style layout.

  The HPS overlay reads one or more local text files whose first two numeric
  columns are interpreted as:
    mass, epsilon

  By default:
    - 2016 HPS files are enabled and drawn in blue.
    - 2021 HPS files are disabled and drawn in red when enabled.
    - 2016 files are searched with exclusion-contour-part-*.txt, excluding
      any filename containing _2021.
    - 2021 files are searched with exclusion-contour-part-*_2021.txt and
      exclusion-contour-part-*.txt_2021.

Inputs:
  data/
    The hps-simpy .dat files copied into this bundle.

  ./
    Optional HPS reach files placed in the same directory as this script, e.g.
      exclusion-contour-part-0.txt
      exclusion-contour-part-1.txt
      exclusion-contour-part-0_2021.txt
      exclusion-contour-part-1_2021.txt

Outputs:
  A PNG, and optionally a PDF, of the styled exclusion plot.

Dependencies:
  Python 3.6+, numpy, scipy, matplotlib.

Quick start:
  python make_paper_style_limits_with_hps_py36.py

Examples:
  # Default behavior: try to draw 2016 files, skip 2021.
  python make_paper_style_limits_with_hps_py36.py

  # Draw both 2016 and 2021.
  python make_paper_style_limits_with_hps_py36.py --plot-2016 --plot-2021

  # Use custom HPS file patterns.
  python make_paper_style_limits_with_hps_py36.py \
    --hps-2016-glob "exclusion-contour-part-*.txt" \
    --hps-2021-glob "exclusion-contour-part-*_2021.txt" \
    --plot-2021

  # If your HPS mass column is in MeV, convert it to GeV.
  python make_paper_style_limits_with_hps_py36.py --hps-mass-scale 0.001

Notes:
  This script does not smooth, interpolate, thin, or downsample the HPS reach
  files. It plots the points in the order they appear in each input file.

Compatibility note:
  This version avoids "from __future__ import annotations", PEP 585 type hints,
  and PEP 604 union syntax so it can run on older Python 3 installations.
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate


# ---------------------------------------------------------------------------
# USER CONFIGURATION
# ---------------------------------------------------------------------------
# These are the two main settings the user asked for.
PLOT_HPS_2016 = True
PLOT_HPS_2021 = False

# Put your HPS text files next to this script, then edit these patterns if the
# names differ from the defaults.
HPS_2016_GLOBS = [
    "exclusion-contour-part-*.txt",
]
HPS_2021_GLOBS = [
    "exclusion-contour-part-*_2021.txt",
    "exclusion-contour-part-*.txt_2021",
]

# The first two numeric columns are used by default.
HPS_MASS_COLUMN = 0
HPS_EPSILON_COLUMN = 1

# HPS contour files from exclusion-contour-part-*.txt use mV[MeV].
# Convert to GeV to match the x-axis m_A' / GeV.
HPS_MASS_SCALE = 0.001

# The HPS contour files may label the second column epsilon^2, but the old
# plot places those numeric values directly on the y-axis. Leave this at
# 1.0 for an exact old-plot overlay. Set to 0.5 only if you want sqrt().
HPS_EPSILON_POWER = 1.0

# HPS overlay styling.
HPS_2016_COLOR = "blue"
HPS_2021_COLOR = "red"
HPS_LINEWIDTH = 2.0
SHOW_HPS_LEGEND = False

# External/reference SIMP plot settings.
DEFAULT_MPIOVERFPI = "4pi"
DEFAULT_LUMI_TEXT = ""
DEFAULT_YEAR_TEXT = "2016"
DEFAULT_OUTPUT = Path("output/hps_style_with_hps.png")


SCAN_FILE_TEMPLATE = "ft_reach_mpi_mrho_ma_1_1.8_3_ad_0.01_mpioverfpi_{mpioverfpi}_ldmx_theory_paper.dat"
FLOAT_RE = re.compile(r"[+-]?(?:(?:\d+\.\d*)|(?:\.\d+)|(?:\d+))(?:[eEdD][+-]?\d+)?")


class PlotInputs(object):
    def __init__(self, data_dir, mpioverfpi):
        self.data_dir = data_dir
        self.mpioverfpi = mpioverfpi

    def path(self, filename):
        return self.data_dir / filename

    def loadtxt(self, filename):
        path = self.path(filename)
        if not path.exists():
            raise FileNotFoundError("Required input file not found: {}".format(path))
        return np.loadtxt(str(path))

    def abundance_files(self):
        if self.mpioverfpi == "4pi":
            return (
                "der_line_rV_1p8_aD_1m2_fine.dat",
                "der_line_rV_1p6_aD_1m2_fine.dat",
            )
        return (
            "abun_cont_mPi_vs_eps_rV_1p8_aD_1m2_neutral_V_only_combined.dat",
            "abun_cont_mPi_vs_eps_rV_1p6_aD_1m2_neutral_V_only.dat",
        )

    def abundances(self):
        first_file, second_file = self.abundance_files()
        first = self.loadtxt(first_file).copy()
        second = self.loadtxt(second_file).copy()

        # hps-simpy/limits.py gives these as m_pi versus epsilon, so convert
        # the x-axis to m_Aprime.
        first[:, 0] = 3.0 * first[:, 0]
        second[:, 0] = 3.0 * second[:, 0]
        return first, second

    def cmb_decay_bound(self):
        return self.loadtxt("cmb_decay_bound_Delta_0p05_mPiOverFpi_{}_aD_1em2.dat".format(self.mpioverfpi))

    def babar(self):
        return self.loadtxt("babar_lim_aD_1em2_mPiOverFpi_{}_two_body.dat".format(self.mpioverfpi))

    def scan_file(self):
        return self.loadtxt(SCAN_FILE_TEMPLATE.format(mpioverfpi=self.mpioverfpi))

    def required_files(self):
        abund_first, abund_second = self.abundance_files()
        return [
            abund_first,
            abund_second,
            "cmb_decay_bound_Delta_0p05_mPiOverFpi_{}_aD_1em2.dat".format(self.mpioverfpi),
            "babar_lim_aD_1em2_mPiOverFpi_{}_two_body.dat".format(self.mpioverfpi),
            SCAN_FILE_TEMPLATE.format(mpioverfpi=self.mpioverfpi),
        ]

    def check_required_files(self):
        missing = [str(self.path(name)) for name in self.required_files() if not self.path(name).exists()]
        if missing:
            missing_lines = "\n".join(missing)
            raise FileNotFoundError("Missing required input files:\n{}".format(missing_lines))


def natural_key(path):
    parts = re.split(r"(\d+)", path.name)
    key = []
    for part in parts:
        if part.isdigit():
            key.append(int(part))
        else:
            key.append(part)
    return key


def log_interpolate_griddata(ref_xy, ref_z, eval_x, eval_y):
    """Replicate the interpolation logic from hps-simpy/limits.py."""
    log_ref_z = np.full_like(ref_z, -200.0)
    positive = ref_z > 0
    log_ref_z[positive] = np.log(ref_z[positive])
    log_grid = interpolate.griddata(
        np.log(ref_xy),
        log_ref_z,
        (np.log(eval_x), np.log(eval_y)),
        method="linear",
    )
    return np.exp(log_grid)


def scan_exclusion_boundary(table, experiment, level, grid_points):
    """Return the single contour boundary used for E137 or Orsay."""
    column_by_experiment = {"e137": 7, "orsay": 8}
    if experiment not in column_by_experiment:
        allowed = ", ".join(column_by_experiment)
        raise ValueError("Unknown experiment '{}'. Allowed values: {}".format(experiment, allowed))

    m_aprime, epsilon = np.meshgrid(
        np.logspace(-2.0, 1.0, grid_points),
        np.logspace(-6.0, -2.0, grid_points),
    )
    interpolated = log_interpolate_griddata(
        table[:, (0, 1)],
        table[:, column_by_experiment[experiment]],
        m_aprime,
        epsilon,
    )

    scratch_fig, scratch_ax = plt.subplots()
    contour = scratch_ax.contour(m_aprime, epsilon, interpolated, levels=[level], linewidths=[0.0])
    if not contour.allsegs or not contour.allsegs[0]:
        plt.close(scratch_fig)
        raise RuntimeError("Could not find contour for {} at level {}".format(experiment, level))

    segments = contour.allsegs[0]
    longest_segment = max(segments, key=len)
    points = np.asarray(longest_segment)[:, (0, 1)]
    plt.close(scratch_fig)
    return points


def fill_scan_exclusion(ax, table, experiment, level, grid_points, y_top, color, zorder):
    boundary = scan_exclusion_boundary(table, experiment, level, grid_points)
    ax.fill_between(
        boundary[:, 0],
        boundary[:, 1],
        y_top,
        color=color,
        linewidth=0.0,
        edgecolor="none",
        zorder=zorder,
    )
    return boundary


def find_hps_files(hps_dir, patterns, year_tag):
    found = []
    for pattern in patterns:
        pattern_path = Path(pattern)
        if pattern_path.is_absolute():
            matches = pattern_path.parent.glob(pattern_path.name)
        elif any(char in pattern for char in "*?[]"):
            matches = hps_dir.glob(pattern)
        else:
            matches = [hps_dir / pattern]

        for path in matches:
            if not path.is_file():
                continue

            # The default 2016 glob also matches names like
            # exclusion-contour-part-0_2021.txt. Keep those out of the
            # 2016 list unless the user gave a more explicit custom pattern.
            if year_tag == "2016" and "_2021" in path.name:
                continue

            found.append(path)

    unique = {}
    for path in found:
        unique[str(path.resolve())] = path
    return sorted(unique.values(), key=natural_key)


def read_hps_reach_file(path, mass_column, epsilon_column, mass_scale, epsilon_power):
    rows = []
    max_column = max(mass_column, epsilon_column)

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.split("#", 1)[0]
            matches = FLOAT_RE.findall(line)
            if not matches:
                continue

            numeric_values = [float(value.replace("D", "E").replace("d", "e")) for value in matches]
            if len(numeric_values) <= max_column:
                continue

            mass = numeric_values[mass_column] * mass_scale
            epsilon = numeric_values[epsilon_column] ** epsilon_power
            rows.append((mass, epsilon))

    if not rows:
        raise ValueError(
            "No usable numeric rows found in {}. Need at least columns {} and {}.".format(
                path,
                mass_column,
                epsilon_column,
            )
        )

    array = np.asarray(rows, dtype=float)
    positive = (array[:, 0] > 0.0) & (array[:, 1] > 0.0)
    if not np.all(positive):
        removed = int(np.size(positive) - np.count_nonzero(positive))
        print("Warning: removed {} non-positive rows from {}; log axes require positive values.".format(removed, path))
        array = array[positive]

    if array.size == 0:
        raise ValueError("All rows in {} were non-positive after filtering for log axes.".format(path))

    return array


def plot_hps_reach_files(
    ax,
    files,
    color,
    label,
    mass_column,
    epsilon_column,
    mass_scale,
    epsilon_power,
    linewidth,
    zorder,
):
    for index, path in enumerate(files):
        points = read_hps_reach_file(
            path=path,
            mass_column=mass_column,
            epsilon_column=epsilon_column,
            mass_scale=mass_scale,
            epsilon_power=epsilon_power,
        )
        line_label = label if index == 0 else None
        ax.plot(
            points[:, 0],
            points[:, 1],
            color=color,
            linewidth=linewidth,
            linestyle="-",
            label=line_label,
            zorder=zorder,
        )
        print(
            "Drew {}: {} ({} points, x=[{:.6g}, {:.6g}], y=[{:.6g}, {:.6g}])".format(
                label,
                path,
                len(points),
                np.min(points[:, 0]),
                np.max(points[:, 0]),
                np.min(points[:, 1]),
                np.max(points[:, 1]),
            )
        )


def add_hps_top_labels(ax, mpioverfpi, lumi_text, year_text):
    title = r"$m_{\pi}/f_{\pi}=4\pi$" if mpioverfpi == "4pi" else r"$m_{\pi}/f_{\pi}=3$"

    ax.text(
        0.0,
        1.025,
        "H P S",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=15,
        fontweight="bold",
    )
    ax.text(
        0.5,
        1.025,
        title,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=15,
    )

    if lumi_text:
        ax.text(
            0.95,
            1.025,
            lumi_text,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=9,
        )

    if year_text:
        ax.text(
            1.015,
            0.98,
            year_text,
            transform=ax.transAxes,
            ha="left",
            va="top",
            rotation=-90.0,
            fontsize=11,
        )

def add_region_labels(ax, mpioverfpi):
    gray_text = "0.45"
    if mpioverfpi == "4pi":
        ax.text(0.42, 2.5e-3, "BaBar", color=gray_text, fontsize=12)
        ax.text(5.0e-2, 1.0e-4, "Orsay", color=gray_text, fontsize=12)
        ax.text(0.13, 2.0e-5, "E137", color=gray_text, fontsize=12)
        ax.text(5.0e-2, 4.5e-6, r"$m_V/m_{\pi}=1.8$", color="black", rotation=18.0, va="bottom", fontsize=12)
        ax.text(5.0e-2, 7.0e-7, r"$m_V/m_{\pi}=1.6$", color="0.35", rotation=18.0, va="bottom", fontsize=12)
        ax.text(1.0, 2.5e-6, "CMB", color=gray_text, fontsize=12)
    else:
        ax.text(0.42, 2.0e-3, "BaBar", color=gray_text, fontsize=12)
        ax.text(5.0e-2, 1.0e-4, "Orsay", color=gray_text, fontsize=12)
        ax.text(0.18, 8.0e-6, "E137", color=gray_text, fontsize=12)
        ax.text(5.0e-2, 3.5e-6, r"$m_V/m_{\pi}=1.8$", color="black", rotation=18.0, va="bottom", fontsize=12)
        ax.text(5.0e-2, 6.0e-7, r"$m_V/m_{\pi}=1.6$", color="0.35", rotation=18.0, va="bottom", fontsize=12)
        ax.text(1.0, 2.5e-6, "CMB", color=gray_text, fontsize=12)

def make_plot(
    data_dir,
    hps_dir,
    output,
    mpioverfpi,
    lumi_text,
    year_text,
    grid_points,
    make_pdf,
    plot_2016,
    plot_2021,
    hps_2016_globs,
    hps_2021_globs,
    hps_mass_column,
    hps_epsilon_column,
    hps_mass_scale,
    hps_epsilon_power,
    hps_2016_color,
    hps_2021_color,
    hps_linewidth,
    show_hps_legend,
):
    plot_inputs = PlotInputs(data_dir=data_dir, mpioverfpi=mpioverfpi)
    plot_inputs.check_required_files()

    x_min = 3.0e-2
    x_max = 2.0
    y_min = 1.0e-7
    y_max = 1.0e-2
    gray_fill = "0.82"

    # Match the old HPS paper-style panel geometry closely.
    fig, ax = plt.subplots(figsize=(7.0, 3.65))

    babar = plot_inputs.babar()
    ax.fill_between(
        babar[:, 0],
        babar[:, 1],
        y_max,
        color=gray_fill,
        linewidth=0.0,
        edgecolor="none",
        zorder=1.0,
    )

    scan = plot_inputs.scan_file()
    fill_scan_exclusion(
        ax=ax,
        table=scan,
        experiment="e137",
        level=10.0,
        grid_points=grid_points,
        y_top=y_max,
        color=gray_fill,
        zorder=1.1,
    )
    fill_scan_exclusion(
        ax=ax,
        table=scan,
        experiment="orsay",
        level=3.0,
        grid_points=grid_points,
        y_top=y_max,
        color=gray_fill,
        zorder=1.2,
    )

    cmb = plot_inputs.cmb_decay_bound()
    ax.plot(cmb[:, 0], cmb[:, 1], color="0.70", linestyle="-", linewidth=1.0, zorder=2.0)

    abundance_18, abundance_16 = plot_inputs.abundances()
    ax.plot(abundance_18[:, 0], abundance_18[:, 1], color="black", linestyle="-", linewidth=2.8, zorder=3.0)
    ax.plot(abundance_16[:, 0], abundance_16[:, 1], color="0.45", linestyle="--", linewidth=2.2, zorder=3.0)

    if plot_2016:
        hps_2016_files = find_hps_files(hps_dir=hps_dir, patterns=hps_2016_globs, year_tag="2016")
        if hps_2016_files:
            plot_hps_reach_files(
                ax=ax,
                files=hps_2016_files,
                color=hps_2016_color,
                label="HPS 2016",
                mass_column=hps_mass_column,
                epsilon_column=hps_epsilon_column,
                mass_scale=hps_mass_scale,
                epsilon_power=hps_epsilon_power,
                linewidth=hps_linewidth,
                zorder=10.0,
            )
        else:
            print(
                "Warning: plot_2016 is enabled, but no 2016 HPS files were found. "
                "Searched in {} with patterns: {}".format(hps_dir, list(hps_2016_globs))
            )

    if plot_2021:
        hps_2021_files = find_hps_files(hps_dir=hps_dir, patterns=hps_2021_globs, year_tag="2021")
        if hps_2021_files:
            plot_hps_reach_files(
                ax=ax,
                files=hps_2021_files,
                color=hps_2021_color,
                label="HPS 2021",
                mass_column=hps_mass_column,
                epsilon_column=hps_epsilon_column,
                mass_scale=hps_mass_scale,
                epsilon_power=hps_epsilon_power,
                linewidth=hps_linewidth,
                zorder=10.1,
            )
        else:
            print(
                "Warning: plot_2021 is enabled, but no 2021 HPS files were found. "
                "Searched in {} with patterns: {}".format(hps_dir, list(hps_2021_globs))
            )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel(r"$m_{A'}$ / GeV", fontsize=16)
    ax.set_ylabel(r"$\epsilon$", fontsize=16)

    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=14, length=8, width=1.2)
    ax.tick_params(which="minor", length=4, width=1.0)
    ax.grid(False)

    add_hps_top_labels(ax=ax, mpioverfpi=mpioverfpi, lumi_text=lumi_text, year_text=year_text)
    add_region_labels(ax=ax, mpioverfpi=mpioverfpi)

    if show_hps_legend:
        ax.legend(loc="lower right", fontsize=7, frameon=False)

    fig.tight_layout(pad=0.45)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output), dpi=300)
    if make_pdf:
        fig.savefig(str(output.with_suffix(".pdf")))
    plt.close(fig)


def add_bool_pair(
    parser,
    positive_name,
    negative_name,
    dest,
    default,
    positive_help,
    negative_help,
):
    group = parser.add_mutually_exclusive_group()
    group.add_argument(positive_name, dest=dest, action="store_true", help=positive_help)
    group.add_argument(negative_name, dest=dest, action="store_false", help=negative_help)
    parser.set_defaults(**{dest: default})


def parse_args():
    parser = argparse.ArgumentParser(description="Make a paper-style HPS/SIMP exclusion plot with optional HPS overlays.")
    parser.add_argument("--data-dir", type=Path, default=Path("data"), help="Directory containing hps-simpy .dat files.")
    parser.add_argument("--hps-dir", type=Path, default=Path("."), help="Directory containing local HPS reach .txt files.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output PNG filename.")
    parser.add_argument("--mpioverfpi", choices=["3", "4pi"], default=DEFAULT_MPIOVERFPI, help="Which m_pi/f_pi data set to draw.")
    parser.add_argument("--lumi-text", default=DEFAULT_LUMI_TEXT, help="Optional top-right luminosity text. Empty string omits it.")
    parser.add_argument("--year-text", default=DEFAULT_YEAR_TEXT, help="Optional rotated year text at the upper-right edge.")
    parser.add_argument("--no-year", action="store_true", help="Omit the rotated year text.")
    parser.add_argument("--grid-points", type=int, default=240, help="Interpolation grid size for E137 and Orsay.")
    parser.add_argument("--pdf", action="store_true", help="Also write a PDF next to the PNG.")

    add_bool_pair(
        parser=parser,
        positive_name="--plot-2016",
        negative_name="--no-plot-2016",
        dest="plot_2016",
        default=PLOT_HPS_2016,
        positive_help="Draw local 2016 HPS reach files.",
        negative_help="Do not draw local 2016 HPS reach files.",
    )
    add_bool_pair(
        parser=parser,
        positive_name="--plot-2021",
        negative_name="--no-plot-2021",
        dest="plot_2021",
        default=PLOT_HPS_2021,
        positive_help="Draw local 2021 HPS reach files.",
        negative_help="Do not draw local 2021 HPS reach files.",
    )

    parser.add_argument(
        "--hps-2016-glob",
        action="append",
        default=None,
        help="Glob or filename for 2016 HPS reach files. Can be repeated.",
    )
    parser.add_argument(
        "--hps-2021-glob",
        action="append",
        default=None,
        help="Glob or filename for 2021 HPS reach files. Can be repeated.",
    )
    parser.add_argument("--hps-mass-column", type=int, default=HPS_MASS_COLUMN, help="Zero-based numeric column index for HPS mass.")
    parser.add_argument("--hps-epsilon-column", type=int, default=HPS_EPSILON_COLUMN, help="Zero-based numeric column index for HPS epsilon.")
    parser.add_argument("--hps-mass-scale", type=float, default=HPS_MASS_SCALE, help="Scale applied to HPS mass before plotting. Default 0.001 converts MeV to GeV.")
    parser.add_argument("--hps-epsilon-power", type=float, default=HPS_EPSILON_POWER, help="Power applied to the HPS epsilon column before plotting. Default 1.0 plots the second column as-is; use 0.5 for sqrt(epsilon^2).")
    parser.add_argument("--hps-2016-color", default=HPS_2016_COLOR, help="Matplotlib color for 2016 HPS contours.")
    parser.add_argument("--hps-2021-color", default=HPS_2021_COLOR, help="Matplotlib color for 2021 HPS contours.")
    parser.add_argument("--hps-linewidth", type=float, default=HPS_LINEWIDTH, help="Line width for HPS contours.")
    parser.add_argument("--hps-legend", action="store_true", default=SHOW_HPS_LEGEND, help="Show a small HPS 2016/2021 legend.")

    return parser.parse_args()


def main():
    args = parse_args()
    year_text = None if args.no_year else args.year_text

    hps_2016_globs = args.hps_2016_glob if args.hps_2016_glob is not None else HPS_2016_GLOBS
    hps_2021_globs = args.hps_2021_glob if args.hps_2021_glob is not None else HPS_2021_GLOBS

    make_plot(
        data_dir=args.data_dir,
        hps_dir=args.hps_dir,
        output=args.output,
        mpioverfpi=args.mpioverfpi,
        lumi_text=args.lumi_text,
        year_text=year_text,
        grid_points=args.grid_points,
        make_pdf=args.pdf,
        plot_2016=args.plot_2016,
        plot_2021=args.plot_2021,
        hps_2016_globs=hps_2016_globs,
        hps_2021_globs=hps_2021_globs,
        hps_mass_column=args.hps_mass_column,
        hps_epsilon_column=args.hps_epsilon_column,
        hps_mass_scale=args.hps_mass_scale,
        hps_epsilon_power=args.hps_epsilon_power,
        hps_2016_color=args.hps_2016_color,
        hps_2021_color=args.hps_2021_color,
        hps_linewidth=args.hps_linewidth,
        show_hps_legend=args.hps_legend,
    )
    print("Wrote {}".format(args.output))
    if args.pdf:
        print("Wrote {}".format(args.output.with_suffix(".pdf")))


if __name__ == "__main__":
    main()

