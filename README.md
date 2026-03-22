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

---

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

    mpic++ -O3 -std=c++11 C++/simulation_2hdm_artificial.cpp -o oscillon_2hdm

---

## 🚀 Quick Start

### 1. Generate oscillon profiles (required for artificial simulations)

    cd Python
    python function_shape.py   # creates 2D_solutions_omega_0.4500_4fields.txt etc.

### 2. Run a domain wall decay simulation (random initial conditions)

    python 2HDM_random.py          # full 2HDM, 2D, random initial fields

To vary the number of active fields, use `4_fields_random.py` or `2_fields_random.py`.  
For 3D, modify the C++ program `mpi_2HDM_random_initial_3D.cpp` and compile.

### 3. Simulate an isolated oscillon (artificial initial condition)

    python evolving-8.py            # 2HDM oscillon with absorbing boundaries

The script reads the previously generated profile file. Output snapshots are saved in `Data/`.

### 4. Perform Floquet stability analysis

    python 2HDM-2D-Floquet.py        # computes eigenvalues and plots Floquet multipliers

For 3D, use `2HDM-3D-Floquet.py` (note: very memory intensive).

### 5. Visualize results
- **2D R1 distribution**: `python 2D_plot.py`
- **GIF animation**: `python plot_gif.py`
- **3D isosurface**: `python 3D_plot.py` (requires 3D data files)

---

## ⏱️ Performance Notes & Running Long Simulations

Some simulations in this repository are computationally demanding. Below are guidelines to help you run them efficiently:

### 🔥 Resource-Intensive Scripts

| Script | Typical Runtime | Memory | Recommended Hardware |
|--------|----------------|--------|----------------------|
| `function_shape.py` | Minutes | < 1 GB | Any desktop |
| `evolving-8.py` (2D, 256² grid) | Hours | 2–4 GB | Desktop / cluster |
| `mpi_2HDM_random_initial_3D.cpp` (128³ grid) | Days | > 32 GB | Cluster with MPI |
| `*Floquet.py` (2D, 129 radial points) | Hours to days | 8–16 GB | High‑memory node |

### 💡 Tips for Reducing Computational Cost

1. **Decrease grid resolution**:  
   - For 2D, change `nx`, `ny` in the script (e.g., from 256 to 128).  
   - For 3D, reducing `nz` will dramatically cut memory/time.

2. **Shorten simulation time**:  
   - Modify `nt` to a smaller value (e.g., half the light crossing time).

3. **Use lower precision in Floquet analysis**:  
   - In `*Floquet.py`, reduce `n_time_samples` (e.g., from 40000 to 5000).  
   - This will still give a qualitative stability picture but with less accuracy.

4. **Run on a cluster / HPC**:  
   - Compile C++ codes with `-O3` optimization and use `mpirun` for parallel runs.  
   - Example for 3D random simulation:  

         mpic++ -O3 -std=c++11 C++/mpi_2HDM_random_initial_3D.cpp -o dw_3d
         mpirun -np 32 ./dw_3d

5. **Use checkpointing**:  
   - The Python scripts do not automatically save intermediate states. For long runs, you may add checkpointing (e.g., save fields every N steps) to resume from a crash.

### 📥 Pre‑computed Data (Optional)

If you only need to reproduce the figures, you can download pre‑computed simulation data from [Zenodo/Figshare link] instead of running the time‑consuming simulations yourself. This data includes:
- Domain wall decay trajectories for various seeds.
- Oscillon profiles for multiple frequencies.
- Floquet eigenvalues and eigenvectors.

<!-- Replace the link with actual DOI/URL when available -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.xxxxxxx.svg)](https://doi.org/10.5281/zenodo.xxxxxxx)

---

## 📊 Reproducing Paper Figures

| Figure | Description | Command / Method |
|--------|-------------|------------------|
| Fig. 1 | Domain wall length evolution | Run `2HDM_random.py` and plot `valsPerLoop*.txt` column 5 |
| Fig. 2 | Lifetime histograms | Run `histogram.py` (requires many seeds) |
| Fig. 3 | Background profiles | `python function_shape.py` |
| Fig. 4 | Oscillon evolution snapshots | `python evolving-8.py` + `python plot_gif.py` |
| Fig. 5 | Floquet multipliers (4‑field) | `python 4-field-2D-Floquet.py` |
| Fig. 6 | Background profiles for 2HDM | `python function_shape.py` (already includes 2HDM) |
| Fig. 7 | Floquet multipliers for 2HDM | `python 2HDM-2D-Floquet.py` |

*For 3D results (Fig. 1b, 3b, etc.), use the corresponding C++ code or the pre‑computed data.*

---

## ⚙️ Parameter Notes

- All parameters (masses, couplings) are set according to Appendices 1 and 2 of the paper.
- Dimensionless units: distances in \(1/m_h\), time in \(1/(c m_h)\) with \(m_h = 125\ \mathrm{GeV}\).
- Random seeds are fixed in Python scripts (`seed=8` for 4‑field, `seed=79` for 2HDM) and can be changed for statistics.
- The C++ codes use hard‑coded parameters; you may modify them at the top of each file.

---

## 🖥️ Running C++ Simulations

For high‑performance runs (especially 3D), compile the C++ code with MPI support:

    mpic++ -O3 -std=c++11 C++/mpi_2HDM_random_initial_2D.cpp -o dw_2hdm
    mpirun -np 4 ./dw_2hdm

**Note**: The MPI versions currently use **1D domain decomposition** (x‑direction). For 3D, a 2D decomposition would be more efficient but is not implemented. Use the serial versions for moderate grid sizes.

---

## 📝 Citation

If you use this code in your research, please cite the accompanying paper:

    @article{meng2025longlived,
      title={Long-Lived Oscillons as Closed Domain Walls in the Z₂-Symmetric Two-Higgs-Doublet Model},
      author={Meng, Zhaoyu},
      journal={arXiv preprint},
      year={2025}
    }

---

## 📬 Contact

For questions or bug reports, please contact **Zhaoyu Meng** (email: [your‑email]) or open an issue on GitHub.

---

## 📜 License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.
