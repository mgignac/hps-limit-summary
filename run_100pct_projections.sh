#!/usr/bin/env bash
set -euo pipefail

# Purpose: Plot the existing HPS exclusion curves (2016, 2021 10%) plus two
# 100%-luminosity projection configurations, styled distinctly so they can
# be compared against each other.
#
# Usage:
#   1. Put your HPS files in the same directory as this script:
#        exclusion-contour-part-0.txt          (2016, dark gray)
#        exclusion-contour-part-1.txt          (2016, dark gray)
#        hps_simps_mu_contour_exp.txt          (2021 10%, solid red)
#        <config A 100% projection file>       (100% proj. A, light-red dashed)
#        <config B 100% projection file>       (100% proj. B, light-blue dashed)
#   2. Set PROJ_A_GLOB and PROJ_B_GLOB below to the filenames for the two
#      100% configurations, and update PROJ_A_LABEL / PROJ_B_LABEL to
#      describe what differs between them.
#   3. Run:
#        bash run_100pct_projections.sh
#
# Notes:
#   - HPS mass columns are in MeV, so --hps-mass-scale 0.001 converts to GeV.
#   - Each projection file must list its points in traced contour order; the
#     plotter connects them in file order and does not reorder.

PROJ_A_GLOB="hps_simps_majd_100_lumi_bkg8.98344x_mu_contour_obs.txt"
PROJ_A_LABEL="HPS 2021 (100% proj., majd)"

PROJ_B_GLOB="hps_simps_rory_100_lumi_bkg8.98344x_mu_contour_obs.txt"
PROJ_B_LABEL="HPS 2021 (100% proj., rory)"

python make_paper_style_limits_with_hps_py36.py \
  --data-dir data \
  --hps-dir . \
  --hps-mass-scale 0.001 \
  --plot-2021 \
  --hps-2021-glob hps_simps_mu_contour_exp.txt \
  --plot-2021-proj \
  --hps-2021-proj-glob "${PROJ_A_GLOB}" \
  --hps-2021-proj-label "${PROJ_A_LABEL}" \
  --plot-2021-proj2 \
  --hps-2021-proj2-glob "${PROJ_B_GLOB}" \
  --hps-2021-proj2-label "${PROJ_B_LABEL}" \
  --hps-legend \
  --no-year \
  --output output/hps_style_with_hps_100pct_projections.png
