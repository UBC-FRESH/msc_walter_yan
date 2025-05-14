# BC Forest Sector Climate Benefits Modeling Framework

This repository contains the mathematical modeling framework developed to simulate and optimize the climate change mitigation potential of the British Columbia (BC) forest sector, encompassing both the forest ecosystem and the harvested wood products (HWP) industry.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Usage](#usage)
- [Data](#data)
- [Utilities](#utilities)
- [License](#license)
- [Contact](#contact)

---

## Overview

The framework combines the [libcbm_py](https://github.com/cat-cfs/libcbm_py) library with a custom [ws3](#) wrapper to:

1. Generate and validate carbon yield curves at stand and regional scales.
2. Optimize harvest and HWP scenarios to minimize net carbon emissions.
3. Perform sensitivity analyses on key HWP parameters.

All analyses are implemented in Jupyter notebooks for reproducibility and extension.

---

## Features

- **Single-Stand Modeling**: Generate, plug in, and validate site-specific carbon yield curves.
- **Regional Modeling**: Scale yield-curve workflows to the full TSA 24 region.
- **Emissions Optimization**: Loop over HWP half-lives and displacement factors to minimize net emissions.
- **Sensitivity Analysis**: Explore how HWP half-life, utilization rate, and displacement factor affect mitigation potential.
- **Modular Utilities**: `util.py` hosts reusable functions for scheduling, coefficient compilation, CBM execution, plotting, and scenario management.

---

## Repository Structure

```text
├── data/
│   ├── woodstock_model_files_test     # SIT and ws3 export for single-stand test
│   └── woodstock_model_files_tsa24    # SIT and ws3 export for TSA 24 region
├── libcbm_py/                         # Local clone of the libcbm_py library
├── ws3/                               # Wrapper modules for ws3 integration
├── generate_c_curves_single_stand.ipynb  # Single-stand carbon curve workflow
├── generate_c_curves_tsa24.ipynb         # TSA 24 regional carbon curve workflow
├── minimize_emission_loop.ipynb          # HWP scenario optimization
├── sensitivity_analysis.ipynb            # Multi-parameter sensitivity analysis
├── util.py                            # Utility functions used by all notebooks
├── LICENSE                            # MIT License
└── README.md                          # This file
