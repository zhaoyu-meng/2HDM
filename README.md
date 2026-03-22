# Long-Lived Oscillons as Closed Domain Walls in the Z₂-Symmetric Two-Higgs-Doublet Model

[![arXiv](https://img.shields.io/badge/arXiv-xxxx.xxxxx-red)](https://arxiv.org/abs/xxxx.xxxxx)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![C++11](https://img.shields.io/badge/C%2B%2B-11-blue.svg)](https://en.cppreference.com/w/cpp/11)

This repository contains the complete source code for the paper:

> **“Long-Lived Oscillons as Closed Domain Walls in the Z₂-Symmetric Two-Higgs-Doublet Model”**  
> Zhaoyu Meng, Huazhong University of Science and Technology, June 2025

The code reproduces all numerical results, including:
- Domain wall decay simulations in 2D and 3D with random initial conditions.
- Reconstruction of radiation‑free oscillon solutions in four‑field and full 2HDM systems.
- Floquet stability analysis of oscillon configurations.
- Generation of publication‑quality figures and animations.

## 📁 Repository Structure
.
├── Data/                          # Output directory (created automatically)
├── Python/                        # Python scripts (main analysis)
│   ├── function_shape.py          # Solve background f0/f1 profiles (4‑field & 2HDM)
│   ├── plot_gif.py                # Generate GIF from simulation snapshots
│   ├── 2D_plot.py                 # Plot 2D R1 distribution and energy density
│   ├── 3D_plot.py                 # Extract and visualize 3D isosurface
│   ├── histogram.py               # Compute and plot bubble lifetimes
│   ├── 4-field-2D-Floquet.py      # Floquet analysis for 4‑field (2D)
│   ├── 4-field-3D-Floquet.py      # Floquet analysis for 4‑field (3D)
│   ├── 2HDM-2D-Floquet.py         # Floquet analysis for 2HDM (2D)
│   └── 2HDM-3D-Floquet.py         # Floquet analysis for 2HDM (3D)
├── C++/                           # High‑performance C++ simulations
│   ├── simulation_2hdm_artificial.cpp
│   ├── simulation_4fields_artificial.cpp
│   ├── 2HDM_random_initial_2D.cpp
│   ├── 2HDM_random_initial_3D.cpp
│   ├── 2HDM_random_initial_loop.cpp
│   ├── 4field_random_initial_loop.cpp
│   └── 2field_random_initial_loop.cpp
├── requirements.txt               # Python dependencies
├── LICENSE                        # MIT License
└── README.md                      # This file
