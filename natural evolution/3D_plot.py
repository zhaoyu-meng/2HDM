import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from skimage import measure
import imageio

# Global parameters
grid_size = 128          # Number of grid points in each dimension
dx = 0.5                 # Grid spacing [1/m_h]
L = grid_size * dx / 2.0 # Half side length of the cubic domain
iso_value = 0.0          # Isosurface value (R1 = 0 corresponds to domain wall)

def load_data(file_name, column_index=1):
    """
    Load data from a file and reshape into a 3D array.

    Parameters
    ----------
    file_name : str
        Path to the data file.
    column_index : int, optional
        Column index (0‑based) to read. Default is 1 (second column), which is R1.

    Returns
    -------
    data_3d : ndarray, shape (grid_size, grid_size, grid_size)
        3D array containing the field values.
    """
    # Read the specified column, skipping the header line
    init_data = np.genfromtxt(file_name, skip_header=1, usecols=column_index)

    # Reshape to a 3D cube
    data_3d = init_data[:grid_size*grid_size*grid_size].reshape(
        (grid_size, grid_size, grid_size)
    )
    return data_3d

def extract_isosurface(data_3d, iso_value=0.0):
    """
    Extract an isosurface from 3D data using the marching cubes algorithm.

    Parameters
    ----------
    data_3d : ndarray
        3D input data.
    iso_value : float, optional
        Isovalue for the surface. Default is 0.0.

    Returns
    -------
    verts, faces, normals, values : tuple
        The output of skimage.measure.marching_cubes.
        Returns (None, None, None, None) if extraction fails.
    """
    try:
        verts, faces, normals, values = measure.marching_cubes(
            data_3d,
            level=iso_value
        )
        return verts, faces, normals, values
    except Exception as e:
        print(f"Error extracting isosurface: {e}")
        return None, None, None, None

def convert_to_physical_coordinates(verts):
    """
    Convert vertex coordinates from grid indices to physical coordinates.

    Grid indices range [0, grid_size-1] and map to physical coordinates
    [-L, L] where L = grid_size * dx / 2.

    Parameters
    ----------
    verts : ndarray
        Vertex coordinates in grid index space.

    Returns
    -------
    verts_physical : ndarray
        Vertex coordinates in physical units.
    """
    if verts is None:
        return None

    verts_physical = np.copy(verts)
    # Shift origin to the center and scale by grid spacing
    for i in range(3):
        verts_physical[:, i] = (verts[:, i] - (grid_size - 1) / 2.0) * dx
    return verts_physical

def create_isosurface_plot(file_name, time_step, column_index=1, save_path=None):
    """
    Create and display/save a 3D isosurface plot of the field.

    Parameters
    ----------
    file_name : str
        Path to the data file.
    time_step : int
        Time step number (used only for title, not saved if save_path is provided).
    column_index : int, optional
        Column index to read from the file (0‑based). Default is 1 (R1).
    save_path : str or None, optional
        If provided, the plot is saved to this path. Otherwise it is displayed.
    """
    # 1. Load the data
    data_3d = load_data(file_name, column_index)
    print(f"Data shape: {data_3d.shape}")
    print(f"Data range: [{data_3d.min():.4f}, {data_3d.max():.4f}]")
    print(f"Physical coordinate range: [-{L:.1f}, {L:.1f}]")

    # 2. Extract the isosurface
    verts, faces, normals, values = extract_isosurface(data_3d, iso_value)
    if verts is None:
        print("Failed to extract isosurface.")
        return
    print(f"Extracted isosurface: {len(verts)} vertices, {len(faces)} faces")

    # 3. Convert to physical coordinates
    verts_physical = convert_to_physical_coordinates(verts)

    # 4. Create the figure
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # 5. Plot the isosurface
    mesh = ax.plot_trisurf(
        verts_physical[:, 0],   # x coordinates (physical)
        verts_physical[:, 1],   # y coordinates (physical)
        faces,                  # triangle connectivity
        verts_physical[:, 2],   # z coordinates (physical)
        cmap='Spectral',
        linewidth=0.5,
        antialiased=True
    )

    # Set axis limits (symmetric about origin)
    ax.set_xlim([-L, L])
    ax.set_ylim([-L, L])
    ax.set_zlim([-L, L])

    # Equal aspect ratio
    ax.set_box_aspect([1, 1, 1])

    # View angle
    ax.view_init(elev=30, azim=45)

    # Grid lines
    ax.grid(True, alpha=0.3)

    # Save or show
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved to: {save_path}")
        plt.close(fig)
    else:
        plt.tight_layout()
        plt.show()

def create_animation(time_steps, base_filename, output_gif="animation.gif"):
    """
    Create a GIF animation of isosurfaces over multiple time steps.

    Parameters
    ----------
    time_steps : list of int
        List of time step numbers to include.
    base_filename : str
        Format string for the file names, e.g.
        "fields_timestep={}_nx128_nt4800_seed8_Z2.txt". The curly braces will be
        replaced by the time step number.
    output_gif : str, optional
        Output GIF file name. Default is "animation.gif".
    """
    frames = []

    for time in time_steps:
        try:
            # Construct the file name
            file_name = base_filename.format(time)

            # Load data and extract isosurface
            data_3d = load_data(file_name, column_index=1)
            verts, faces, normals, values = extract_isosurface(data_3d, iso_value)

            if verts is None:
                continue

            # Convert to physical coordinates
            verts_physical = convert_to_physical_coordinates(verts)

            # Create the plot
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')

            mesh = ax.plot_trisurf(
                verts_physical[:, 0],
                verts_physical[:, 1],
                faces,
                verts_physical[:, 2],
                cmap='Spectral',
                linewidth=0.2,
                antialiased=True
            )

            # Set labels and title
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_title(f'Time Step: {time}')

            # Set axis limits
            ax.set_xlim([-L, L])
            ax.set_ylim([-L, L])
            ax.set_zlim([-L, L])

            ax.set_box_aspect([1, 1, 1])
            ax.view_init(elev=30, azim=45)

            plt.tight_layout()

            # Convert the figure to an image array
            fig.canvas.draw()
            image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
            image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
            frames.append(image)

            plt.close(fig)

            print(f"Processed time step {time}")

        except Exception as e:
            print(f"Error processing time step {time}: {e}")
            continue

    # Save the GIF
    if frames:
        imageio.mimsave(output_gif, frames, fps=5)
        print(f"Animation saved to: {output_gif}")
    else:
        print("No frames available to create animation.")

# Main program
if __name__ == "__main__":
    # Example: generate a single isosurface plot
    print("=== Generating a single time step isosurface ===")
    print(f"Grid size: {grid_size} × {grid_size} × {grid_size}")
    print(f"Grid spacing dx = {dx}")
    print(f"Physical coordinate range: [-{L:.2f}, {L:.2f}]")

    file_name = 'Data/fields_timestep=3600_3d_2hdm_nx_0_128_nt4800_random.txt'

    # Ensure the file exists before trying to plot
    try:
        create_isosurface_plot(
            file_name=file_name,
            time_step=0,
            column_index=1,
            save_path="3D_isosurface.png"   # Set to None to display instead
        )
    except FileNotFoundError:
        print(f"File not found: {file_name}")
        print("Please make sure the file path is correct or use example data.")

    # Uncomment the following block to create an animation
    """
    print("\n=== Creating animation ===")
    time_steps = [0, 400, 800, 1200, 1600, 2000]
    base_filename = "fields_timestep={}_nx128_nt4800_seed8_Z2.txt"
    try:
        create_animation(time_steps, base_filename, "isosurface_animation.gif")
    except Exception as e:
        print(f"Error creating animation: {e}")
    """
