# HPS/SIMP paper-style limit plot with optional local HPS contours

This is the Python 3.6-compatible version. It removes `from __future__ import annotations` and newer type-hint syntax.

## Default behavior

```python
PLOT_HPS_2016 = True
PLOT_HPS_2021 = False
```

## Expected local HPS files

Copy your HPS text files into this same directory after unzipping.

Default 2016 pattern:

```text
exclusion-contour-part-*.txt
```

Default 2021 patterns:

```text
exclusion-contour-part-*_2021.txt
exclusion-contour-part-*.txt_2021
```

The script uses the first two numeric columns as `(mass, epsilon)` and plots the points in the order they appear in each file. It does not smooth, interpolate, thin, or downsample the HPS contours.

## Run

```bash
python -m pip install -r requirements.txt
python make_paper_style_limits_with_hps_py36.py
```

To plot both 2016 and 2021:

```bash
python make_paper_style_limits_with_hps_py36.py --plot-2016 --plot-2021
```

If the HPS mass is in MeV rather than GeV:

```bash
python make_paper_style_limits_with_hps_py36.py --hps-mass-scale 0.001
```
