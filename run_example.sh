#!/usr/bin/env bash
set -euo pipefail

# README:
# Purpose: Run the default HPS/SIMP exclusion plot.
# Usage:
#   1. Unzip this bundle.
#   2. Copy your HPS files into the same directory as this script:
#        exclusion-contour-part-0.txt          (2016, blue)
#        exclusion-contour-part-1.txt          (2016, blue)
#        hps_simps_mu_contour_exp.txt          (2021 10%, solid red)
#        hps_simps_100_lumi_obs_rate_sig8.98339x_obsrate8.98339x_bkgsyst_mu_contour_obs.txt
#                                              (2021 100% projection, light-red dotted)
#   3. Run:
#        bash run_example.sh
#
# Notes:
#   - HPS mass columns are in MeV, so --hps-mass-scale 0.001 converts to GeV.
#   - The 2021 file must list its points in traced contour order; the plotter
#     connects them in file order and does not reorder.

python make_paper_style_limits_with_hps_py36.py \
  --data-dir data \
  --hps-dir . \
  --hps-mass-scale 0.001 \
  --plot-2021 \
  --hps-2021-glob hps_simps_4pi_100_lumi_bkg8.98344x_mu_contour_obs.txt \
  --plot-2021-proj2 \
  --hps-2021-proj2-glob hps_simps_4pi_100_lumi_bkg8.98344x_mu_contour_obs_bugfix.txt \
  --hps-legend \
  --no-year \
  --output output/hps_style_with_hps.png


#   --hps-2021-proj-glob hps_simps_100_lumi_obs_rate_sig8.98339x_obsrate8.98339x_bkgsyst_mu_contour_obs.txt \
