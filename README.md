# BC Forest Sector Climate Benefits Modeling Framework

This repository contains the mathematical modeling framework developed to simulate and optimize the climate change mitigation potential of the British Columbia (BC) forest sector, encompassing both the forest ecosystem and the harvested wood products (HWP) industry.

---

## Overview

The framework integrates [libcbm_py](https://github.com/cat-cfs/libcbm_py) with [ws3](#)  to:

1. Generate, plug-in, validate carbon yield curve method at stand and lanscape levels.
2. Compare and find optimize harvest and HWP scenario to maximize the total carbon stock and minimize net carbon emissions from the forest system.
3. Perform sensitivity analyses on the impact of key HWP climate parameters on the optimal forest carbon management plans.

All analyses are implemented in Jupyter notebooks for reproducibility and extension.


---

## Repository Structure

```text
├── data/
│   ├── woodstock_model_files_test     # Files to build singe-stand test ws3 model
│   └── woodstock_model_files_tsa24    # Files to build TSA 24 ws3 model
├── libcbm_py/                         # Local clone of the libcbm_py library
├── ws3/                               # Local clone of the ws3 library
├── generate_c_curves_single_stand.ipynb  # Single-stand carbon yield curve generation, plug-in, and validation workflow
├── generate_c_curves_tsa24.ipynb         # TSA 24 lanscape carbon yield curve generation, plug-in, and validation workflow
├── minimize_emission_loop.ipynb          # Loop to model each minimize-system-emission scenario with increasing HWP half-life and displacement factor
├── sensitivity_analysis.ipynb            # Analyze the sensitivity of optimal modelling solution to HWP half-life and displacement factor
├── util.py                            # Utility functions used by all notebooks
├── LICENSE                            # MIT License
└── README.md                          # This file
