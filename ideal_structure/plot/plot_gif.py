"""
plot_gif.py
------------
Generate GIF animations and field snapshots from simulation data.
This script reads the field snapshot files produced by the C++ simulation
and creates publication‑quality plots of the R1 distribution, radial field
profiles, and energy density. The output can be saved as a GIF animation.

Author: Zhaoyu Meng
Date:   
License: MIT
"""

import numpy as np
import matplotlib.pyplot as plt
import imageio
import os
import io

# Simulation parameters (must match the C++ simulation)
size = 256       # grid points per dimension (nx = ny)
dx = 0.25        # grid spacing in physical units [1/m_h]
dt = 0.1         # time step [1/(c m_h)]

# ----------------------------------------------------------------------
def finalplot(field_indices, field_labels, time):
    """
    Create a figure with three subplots:
      1. 2D color map of R1 (the order parameter used to identify domain walls)
      2. Radial profile of selected fields along the horizontal midline
      3. 2D color map of the energy density

    Parameters
    ----------
    field_indices : list of int
        Column indices (0‑based) of the fields to plot in the radial profile.
    field_labels : list of str
        LaTeX labels for the fields (e.g., "$\\psi_1$").
    time : int
        Time step index (the snapshot file name contains this number).

    Returns
    -------
    fig : matplotlib.figure.Figure or None
        The generated figure, or None if the file cannot be read.
    """
    plot_size = 64          # number of grid points to show around the centre
    half = size // 2

    # Create figure with three subplots
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    ax1, ax2, ax3 = axes

    # Build file name: Data/fields_timestep=<time>_2d_2hdm_nx_0_256_nt320_stable.txt
    filename = f'Data/fields_timestep={time}_2d_2hdm_nx_0_256_nt320_stable.txt'
    try:
        data = np.genfromtxt(filename, skip_header=1)
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        return None

    # Physical coordinates for the plot window
    x_min, x_max = -plot_size * dx, plot_size * dx
    y_min, y_max = -plot_size * dx, plot_size * dx

    # --- Subplot 1: R1 distribution ---
    # R1 is stored in column 1 (0‑based index 1)
    r1_field = data[:size * size, 1].reshape((size, size))
    r1_crop = r1_field[half - plot_size:half + plot_size,
                       half - plot_size:half + plot_size]
    im1 = ax1.imshow(r1_crop, cmap='seismic', aspect='auto',
                     vmin=-2.5, vmax=2.5,
                     extent=[x_min, x_max, y_min, y_max])
    plt.colorbar(im1, ax=ax1)
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_title(f'R1 distribution at t={time*dt:.2f}')

    # --- Subplot 2: Radial field profiles ---
    # Take a horizontal cut through the centre (y = 0)
    r_values = np.linspace(0, plot_size * dx, plot_size)
    for idx, col in enumerate(field_indices):
        label = field_labels[idx]
        if col < data.shape[1]:
            field_data = data[:size * size, col].reshape((size, size))
            # Extract horizontal line from centre to right edge
            profile = field_data[half, half:half + plot_size]
            ax2.plot(r_values, profile, linewidth=1, label=label)

    ax2.set_xlabel("r")
    ax2.set_ylabel("Field value")
    ax2.set_ylim(-2, 2)
    ax2.legend(loc='upper right', fontsize='small')
    ax2.set_title('Field values')

    # --- Subplot 3: Energy density ---
    # Energy density is stored in the last column (index -3)
    energy = data[:size * size, -3].reshape((size, size))
    energy_crop = energy[half - plot_size:half + plot_size,
                         half - plot_size:half + plot_size]
    im3 = ax3.imshow(energy_crop, cmap='seismic', aspect='auto',
                     vmin=-1, vmax=1,
                     extent=[x_min, x_max, y_min, y_max])
    plt.colorbar(im3, ax=ax3)
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")
    ax3.set_title(f'Energy density at t={time*dt:.2f}')

    plt.tight_layout()
    return fig


# ----------------------------------------------------------------------
def create_gif(field_indices, field_labels, start_time, end_time,
               time_step=10, gif_name='animation.gif'):
    """
    Create a GIF animation by looping over time steps.

    Parameters
    ----------
    field_indices, field_labels : see finalplot()
    start_time, end_time : int
        First and last time step indices (inclusive).
    time_step : int
        Interval between frames.
    gif_name : str
        Output file name for the GIF.
    """
    images = []
    print(f"Creating GIF from time={start_time} to {end_time} with step={time_step}")

    # Create output directory if it doesn't exist
    if not os.path.exists('gif'):
        os.makedirs('gif')

    for t in range(start_time, end_time + 1, time_step):
        print(f"Processing time step: {t}")
        fig = finalplot(field_indices, field_labels, t)
        if fig is not None:
            # Save figure to an in‑memory PNG buffer
            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=100)
            buf.seek(0)
            image = imageio.imread(buf)
            images.append(image)
            plt.close(fig)       # free memory
            buf.close()

    if images:
        print(f"Saving GIF with {len(images)} frames...")
        imageio.mimsave(gif_name, images, fps=20, loop=0)
        print(f"GIF saved as {gif_name}")
    else:
        print("No images to save.")


# ----------------------------------------------------------------------
# Main script
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # List of field indices and their LaTeX labels (8 real fields)
    field_indices = [9, 10, 13, 14, 11, 12, 15, 16]
    field_labels = ["$\\psi_1$", "$\\psi_2$", "$\\psi_3$", "$\\psi_4$",
                    "$\\psi_5$", "$\\psi_6$", "$\\psi_7$", "$\\psi_8$"]
    
    # may swap to plot 4-field system by using parameter:
    # field_indices = [9, 10, 13, 14]
    # field_labels = ["$\psi_1$", "$\psi_2$", "$\psi_3$", "$\psi_4$"]

    # Plot a few single frames for inspection
    for t in [0, 35, 70, 105]:
        fig = finalplot(field_indices, field_labels, t)
        if fig:
            plt.show()
            plt.close(fig)

    # Generate the full animation (adjust time range as needed)
    create_gif(field_indices, field_labels,
               start_time=0, end_time=315, time_step=5,
               gif_name='field_evolution.gif')
