# Long-Lived Oscillons as Closed Domain Walls in the Z₂-Symmetric Two-Higgs-Doublet Model

[![arXiv](https://img.shields.io/badge/arXiv-2601.14676-red)](https://arxiv.org/abs/2601.14676)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![C++11](https://img.shields.io/badge/C%2B%2B-11-blue.svg)](https://en.cppreference.com/w/cpp/11)

This repository contains the complete source code for the paper:

> **“Long-Lived Oscillons as Closed Domain Walls in the Z₂-Symmetric Two-Higgs-Doublet Model”**  
> Zhaoyu Meng, Mar 2026
> arXiv:2601.14676 [hep-ph]
The code reproduces all numerical results, including:
- Domain wall decay simulations in 2D and 3D with random initial conditions.
- Reconstruction of radiation‑free oscillon solutions in four‑field and full 2HDM systems.
- Floquet stability analysis of oscillon configurations.
- Generation of publication‑quality figures and animations.

---

## 📁 Repository Structure

    .
    ├── Data/                          # Output directory (created automatically)
    ├── natural_evolution/             # Analysis for random initial condition
    │   ├── 2HDM_random_initial_2D.cpp
    │   ├── 2HDM_random_initial_3D.cpp
    │   ├── 2HDM_random_initial_loop.cpp
    │   ├── 4field_random_initial_loop.cpp
    │   ├── 2field_random_initial_loop.cpp
    │   ├── 2D_plot.py                 # Plot 2D R1 distribution
    │   ├── 3D_plot.py                 # Extract and visualize 3D isosurface
    │   ├── histogram.py               # Compute and plot bubble lifetimes
    ├── ideal_structure/               # Reproduction of the long-lived structure
    │   ├──plot/
    │   |   ├── simulation_2hdm_artificial.cpp
    │   |   ├── simulation_4fields_artificial.cpp
    │   |   ├── function_shape.py          # Solve background f0/f1 profiles (4‑field & 2HDM)
    │   |   └── plot_gif.py                # Generate GIF from simulation snapshots
    │   └── Floquet/
    │       ├── 2HDM-2D-Floquet.py         # Floquet analysis for 2HDM (2D)
    │       ├── 2HDM-3D-Floquet.py         # Floquet analysis for 2HDM (3D)
    │       ├── 4-field-2D-Floquet.py      # Floquet analysis for 4‑field (2D)
    │       └── 4-field-3D-Floquet.py      # Floquet analysis for 4‑field (3D)   
    ├── requirements.txt               # Python dependencies
    ├── LICENSE                        # MIT License
    ├── 2D_solutions_omega_0.45.txt    # output data of f0/f1 values, used for artificial setup
    └── README.md                      # This file

---

## 📦 Dependencies

### Python
- Python 3.8+
- numpy, scipy, matplotlib, imageio, scikit-image

Install with:

    pip install -r requirements.txt

### C++
- C++11 compiler (g++, icc, etc.)
- MPI library (OpenMPI/MPICH) for parallel runs

Compile example:

    cd natural_evolution
    mpic++ -O3 -std=c++11 2HDM_random_initial_2D.cpp -o dw_2hdm

---
## 🚀 Quick Start

### Natural evolution

#### 1. Run a domain wall decay simulation (random initial conditions)

    cd natural_evolution
    ./2HDM_random_initial_2D          # 2D, full 2HDM, start with random field distribution

### Reproduction of the long-lived structure

#### 1. Generate background oscillon profiles

    cd ideal_structure/plot
    python function_shape.py

This creates `2D_solutions_omega_0.45.txt` in the root directory.

#### 2. Run an artificial oscillon simulation

    cd ideal_structure/plot
    ./simulation_2hdm_artificial      # 2HDM oscillon with absorbing boundaries

For the 4‑field system, compile and run `simulation_4fields_artificial.cpp`.  
The script reads `2D_solutions_omega_0.45.txt`. Output snapshots are saved in `Data/`.

#### 3. Perform Floquet stability analysis

    cd ideal_structure/Floquet
    python 2HDM-2D-Floquet.py         # computes eigenvalues and plots Floquet multipliers

For 3D, use `2HDM-3D-Floquet.py`.

#### 4. Visualize results

- **GIF animation**: `cd ideal_structure/plot && python plot_gif.py`
---

## ⏱️ Performance Notes & Running Long Simulations

Some simulations in this repository are computationally demanding. Below are guidelines to help you run them efficiently:

### 🔥 Resource-Intensive Scripts

| Script | Typical Runtime | Memory | Recommended Hardware |
|--------|----------------|--------|----------------------|
| `function_shape.py` | Minutes | < 1 GB | Any desktop |
| `evolving-8.py` (2D, 256² grid) | Hours | 2–4 GB | Desktop / cluster |
| `mpi_2HDM_random_initial_3D.cpp` (128³ grid) | Days | 8–16 GB | High‑memory node  |
| `*Floquet.py` (2D, 129 radial points) | Hours | 2–4 GB | Desktop / cluster |

### 💡 Tips for Reducing Computational Cost
but with less accuracy

1. **Decrease grid resolution**:  
   - For 2D, change `nx`, `ny` in the script (e.g., from 256 to 128).  

2. **Shorten simulation time**:  
   - Modify `nt` to a smaller value.

3. **Use lower precision in Floquet analysis**:  
   - In `*Floquet.py`, reduce `n_time_samples` (e.g., from 40000 to 5000) or number of radial points `Nr`.  
   - This will still give a qualitative stability picture but with less accuracy.

4. **Run on a cluster / HPC**:  
   - Compile C++ codes with `-O3` optimization and use `mpirun` for parallel runs.  
   - Example for 3D random simulation:  

         mpic++ -O3 -std=c++11 C++/mpi_2HDM_random_initial_3D.cpp -o dw_3d
         mpirun -np 32 ./dw_3d

---

## 📊 Reproducing Paper Figures

| Figure | Description | Command / Method |
|--------|-------------|------------------|
| Fig. 1 | Domain wall length evolution | Run `2HDM_random.py` and plot `valsPerLoop*.txt` column 5 |
| Fig. 2 | Lifetime histograms | Run `2HDM_random_initial_loop.cpp` + `4fields_random_initial_loop.cpp` + `2fields_random_initial_loop.cpp` then `histogram.py` (requires many seeds) |
| Fig. 3 | Background profiles | `python function_shape.py` |
| Fig. 4 | Oscillon evolution snapshots | `./simulation_4fields_artificial` + `python plot_gif.py` |
| Fig. 5 | Floquet multipliers (4‑field) | `python 4-field-2D-Floquet.py` +  `python 4-field-3D-Floquet.py` |
| Fig. 6 | Floquet multipliers for 2HDM | `python 2HDM-2D-Floquet.py` + `python 2HDM-3D-Floquet.py` |

*For 3D results (Fig. 1b, 3b, etc.), use the corresponding C++ code or the pre‑computed data.*

---

## ⚙️ Parameter Notes

- All parameters (masses, couplings) are set according to Appendices 1 and 2 of the paper.
- Dimensionless units: distances in \(1/m_h\), time in \(1/(c m_h)\) with \(m_h = 125\ \mathrm{GeV}\).
- Random seeds are fixed in Python scripts (`seed=8` for 3D 2HDM, `seed=79` for 2D 2HDM) and can be changed for statistics.
- The C++ codes use hard‑coded parameters; you may modify them at the top of each file.

---

## 🖥️ Running C++ Simulations

For high‑performance runs (especially 3D), compile the C++ code with MPI support:

    cd natural_evolution
    mpic++ -O3 -std=c++11 2HDM_random_initial_2D.cpp -o dw_2hdm
    mpirun -np 4 ./dw_2hdm

**Note**: The MPI versions currently use **1D domain decomposition** (x‑direction). For 3D, a 2D decomposition would be more efficient but is not implemented. Use the serial versions for moderate grid sizes.

---

## 📝 Citation

If you use this code in your research, please cite the accompanying paper:

    @misc{meng2026longlivedoscillonscloseddomain,
      title={Long-Lived Oscillons as Closed Domain Walls in the $\mathbb Z_2$-Symmetric Two-Higgs-Doublet Model}, 
      author={Zhaoyu Meng},
      year={2026},
      eprint={2601.14676},
      archivePrefix={arXiv},
      primaryClass={hep-ph},
      url={https://arxiv.org/abs/2601.14676}, 
    }

---

## 📬 Contact

For questions or bug reports, please contact **Zhaoyu Meng** (email: [your‑email]) or open an issue on GitHub.

---

## 📜 License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.
