#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
function_shape.py

Solve the background equations for the oscillon profiles f₀(r) and f₁(r)
in the 4‑field model (two complex scalars) under the radiation‑free
parameter choice. The equations correspond to Eq. (4.5) of the paper.

"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid


# ----------------------------------------------------------------------
#  Solver for 2D polar coordinates
# ----------------------------------------------------------------------
def solve_2d_polar(omega, f0_initial, f1_initial, R=32.0, N=1000,
                   mu=0.5, lamb=1.78, lamb4=-2.56):
    """
    Solve the coupled ODEs for f₀(r) and f₁(r) in 2D polar coordinates.

    Parameters
    ----------
    omega : float
        Oscillon frequency (in units of m_h).
    f0_initial, f1_initial : array_like
        Initial guesses for f₀(r) and f₁(r) on the radial grid.
    R : float
        Maximum radius for the solution.
    N : int
        Number of radial grid points.
    mu, lamb, lamb4 : float
        Dimensionless potential parameters (see paper Appendix 1).

    Returns
    -------
    r : ndarray
        Radial grid points (0 to R).
    f0_solution, f1_solution : ndarray
        Solved profiles f₀(r) and f₁(r).
    """
    dr = R / N
    r = np.linspace(0, R, N + 1)
    r_safe = r.copy()
    r_safe[0] = 1e-11

    def equations(y):
        f0 = y[:N+1]
        f1 = y[N+1:2*(N+1)]
        res0 = np.zeros(N+1)
        res1 = np.zeros(N+1)

        # f0 equation (2D polar Laplacian)
        for i in range(N+1):
            if i == 0:
                laplacian = 2 * (f0[1] - f0[0]) / dr**2
            elif i == N:
                laplacian = (2 * f0[N-1] - 2 * f0[N]) / dr**2
            else:
                second = (f0[i-1] - 2*f0[i] + f0[i+1]) / dr**2
                first  = (f0[i+1] - f0[i-1]) / (2 * dr)
                laplacian = second + (1 / r_safe[i]) * first
            nonlinear = mu * f0[i] - 2*lamb * f0[i] * (2*f0[i]**2 + 2*f1[i]**2) \
                        + 2*lamb4 * f0[i] * (-f0[i]**2 + f1[i]**2)
            res0[i] = laplacian + nonlinear

        # f1 equation
        for i in range(N+1):
            if i == 0:
                laplacian = 2 * (f1[1] - f1[0]) / dr**2
            elif i == N:
                laplacian = (2 * f1[N-1] - 2 * f1[N]) / dr**2
            else:
                second = (f1[i-1] - 2*f1[i] + f1[i+1]) / dr**2
                first  = (f1[i+1] - f1[i-1]) / (2 * dr)
                laplacian = second + (1 / r_safe[i]) * first
            nonlinear = omega**2 * f1[i] + mu * f1[i] \
                        - 2*lamb * f1[i] * (2*f0[i]**2 + 2*f1[i]**2) \
                        - 2*lamb4 * f1[i] * (-f0[i]**2 + f1[i]**2)
            res1[i] = laplacian + nonlinear

        return np.concatenate([res0, res1])

    y_init = np.concatenate([f0_initial, f1_initial])
    print(f"Solving 2D system (ω={omega})...")
    sol = fsolve(equations, y_init, xtol=1e-13, maxfev=10_000_000)
    print("2D solution found.")
    return r, sol[:N+1], sol[N+1:]


# ----------------------------------------------------------------------
#  Solver for 3D spherical coordinates
# ----------------------------------------------------------------------
def solve_3d_spherical(omega, f0_initial, f1_initial, R=32.0, N=1000,
                       mu=0.5, lamb=1.78, lamb4=-2.56):
    """
    Solve the coupled ODEs for f₀(r) and f₁(r) in 3D spherical coordinates.

    The Laplacian in 3D is ∇²f = f'' + (2/r) f'.
    """
    dr = R / N
    r = np.linspace(0, R, N + 1)
    r_safe = r.copy()
    r_safe[0] = 1e-11

    def equations(y):
        f0 = y[:N+1]
        f1 = y[N+1:2*(N+1)]
        res0 = np.zeros(N+1)
        res1 = np.zeros(N+1)

        for i in range(N+1):
            if i == 0:
                laplacian0 = 3 * (f0[1] - f0[0]) / dr**2
                laplacian1 = 3 * (f1[1] - f1[0]) / dr**2
            elif i == N:
                laplacian0 = (2 * f0[N-1] - 2 * f0[N]) / dr**2
                laplacian1 = (2 * f1[N-1] - 2 * f1[N]) / dr**2
            else:
                second0 = (f0[i-1] - 2*f0[i] + f0[i+1]) / dr**2
                first0  = (f0[i+1] - f0[i-1]) / (2 * dr)
                laplacian0 = second0 + (2 / r_safe[i]) * first0

                second1 = (f1[i-1] - 2*f1[i] + f1[i+1]) / dr**2
                first1  = (f1[i+1] - f1[i-1]) / (2 * dr)
                laplacian1 = second1 + (2 / r_safe[i]) * first1

            nonlinear0 = mu * f0[i] - 2*lamb * f0[i] * (2*f0[i]**2 + 2*f1[i]**2) \
                        + 2*lamb4 * f0[i] * (-f0[i]**2 + f1[i]**2)
            nonlinear1 = omega**2 * f1[i] + mu * f1[i] \
                        - 2*lamb * f1[i] * (2*f0[i]**2 + 2*f1[i]**2) \
                        - 2*lamb4 * f1[i] * (-f0[i]**2 + f1[i]**2)

            res0[i] = laplacian0 + nonlinear0
            res1[i] = laplacian1 + nonlinear1

        return np.concatenate([res0, res1])

    y_init = np.concatenate([f0_initial, f1_initial])
    print(f"Solving 3D system (ω={omega})...")
    sol = fsolve(equations, y_init, xtol=1e-13, maxfev=10_000_000)
    print("3D solution found.")
    return r, sol[:N+1], sol[N+1:]


# ----------------------------------------------------------------------
#  Utility functions for derivatives and integrals
# ----------------------------------------------------------------------
def derivative(f, r, dr):
    """Compute df/dr using central differences, handling boundaries."""
    N = len(r) - 1
    df = np.zeros_like(r)
    for i in range(N+1):
        if i == 0:
            df[i] = (f[1] - f[0]) / dr
        elif i == N:
            df[i] = (f[N] - f[N-1]) / dr
        else:
            df[i] = (f[i+1] - f[i-1]) / (2 * dr)
    return df


def compute_Q(r, f1, omega, dim):
    """
    Compute Q = ∫ 2 ω f₁² dV.
    dim = 2 for polar, 3 for spherical.
    """
    integrand = 2 * omega * f1**2
    if dim == 2:
        dV = 2 * np.pi * r
    elif dim == 3:
        dV = 4 * np.pi * r**2
    else:
        raise ValueError("dimension must be 2 or 3")
    return trapezoid(integrand * dV, r)


def compute_E(r, f0, f1, omega, mu, lamb, lamb4, dim):
    """
    Compute total energy E from Eq. (4.5) integrated over space.
    """
    dr = r[1] - r[0]
    df0 = derivative(f0, r, dr)
    df1 = derivative(f1, r, dr)

    term1 = (omega**2 * f1)**2
    term2 = df0**2
    term3 = df1**2
    term4 = -mu * (f0**2 + f1**2)
    term5 = 2 * lamb * (f0**2 + f1**2)**2
    term6 = lamb4 * (-f0**2 + f1**2)**2

    integrand = (term1 + term2 + term3 + term4 + term5 + term6 + 0.0625) * 2

    if dim == 2:
        dV = 2 * np.pi * r
    else:
        dV = 4 * np.pi * r**2
    return trapezoid(integrand * dV, r)


# ----------------------------------------------------------------------
#  Main script: compute profiles for a set of frequencies (including ω = 0.45)
# ----------------------------------------------------------------------
def main():
    # Parameters (same as in the paper)
    mu = 0.5
    lamb = 1.78
    lamb4 = -2.56
    R = 32.0
    N = 320
    r = np.linspace(0, R, N+1)

    # List of frequencies and corresponding initial guesses (2D and 3D)
    # (The guesses are tuned to converge; many are commented out to keep the
    #  script focused on the main frequency ω = 0.45 used in the paper.)
    omega_configs = [
        (0.3,
         ((1+np.tanh((r-10.0)/1.5))*0.25, (1-np.tanh((r-9.5)/1.5))*0.27),
         ((1+np.tanh((r-19.0)/1.5))*0.25, (1-np.tanh((r-18.5)/1.5))*0.26)),

        (0.37,
         ((1+np.tanh((r-6.5)/2))*0.25, (1-np.tanh((r-6.5)/2))*0.27),
         ((1+np.tanh((r-12.5)/2))*0.25, (1-np.tanh((r-12.5)/2))*0.26)),

        (0.45,
         ((1+np.tanh((r-5.0)/1.5))*0.25, (1-np.tanh((r-4.5)/1.5))*0.29),
         ((1+np.tanh((r-9.0)/1.5))*0.25, (1-np.tanh((r-8.5)/1.5))*0.29)),

        (0.6,
         ((1+np.tanh((r-3.5)/2))*0.25, (1-np.tanh((r-3.5)/2))*0.29),
         ((1+np.tanh((r-5)/2))*0.25, (1-np.tanh((r-5)/2))*0.31)),

        (0.75,
         ((1+np.tanh((r-3)/2))*0.25, (1-np.tanh((r-3)/2))*0.29),
         ((1+np.tanh((r-3)/2))*0.25, (1-np.tanh((r-3)/2))*0.29)),

        (1.0,
         (0.5 * np.ones(N+1) - 0.4 * (np.cosh(r/1.5)**(-1)), 0.4 * (np.cosh(r/1.5)**(-1))),
         (0.5 * np.ones(N+1) - 0.5 * (np.cosh(r/1.5)**(-1)), 0.6 * (np.cosh(r/1.5)**(-1)))),

        (1.2,
         (0.5 * np.ones(N+1) - 0.3 * (np.cosh(r/0.5)**(-1)), 0.35 * (np.cosh(r/0.5)**(-1))),
         (0.5 * np.ones(N+1) - 0.35 * (np.cosh(r/1)**(-1)), 0.5 * (np.cosh(r/1)**(-1)))),

        (1.4,
         (0.5 * np.ones(N+1) - 0.15 * (np.cosh(r/1.5)**(-1)), 0.2 * (np.cosh(r/1.5)**(-1))),
         (0.5 * np.ones(N+1) - 0.25 * (np.cosh(r/1)**(-1)), 0.35 * (np.cosh(r/1)**(-1)))),

        (1.5,
         (0.5 * np.ones(N+1) - 0.1 * (np.cosh(r/2.5)**(-1)), 0.1 * (np.cosh(r/2.5)**(-1))),
         (0.5 * np.ones(N+1) - 0.1 * (np.cosh(r/2.5)**(-1)), 0.1 * (np.cosh(r/2.5)**(-1)))),

        (1.53,
         (0.5 * np.ones(N+1) - 0.1 * (np.cosh(r/2.5)**(-1)), 0.1 * (np.cosh(r/2.5)**(-1))),
         (0.5 * np.ones(N+1) - 0.1 * (np.cosh(r/2.5)**(-1)), 0.1 * (np.cosh(r/2.5)**(-1)))),

        (1.55,
         (0.5 * np.ones(N+1) - 0.1 * (np.cosh(r/3)**(-1)), 0.1 * (np.cosh(r/3)**(-1))),
         (0.5 * np.ones(N+1) - 0.1 * (np.cosh(r/3)**(-1)), 0.1 * (np.cosh(r/3)**(-1)))),

        (1.57,
         (0.5 * np.ones(N+1) - 0.04 * (np.cosh(r/5)**(-1)), 0.05 * (np.cosh(r/5)**(-1))),
         (0.5 * np.ones(N+1) - 0.04 * (np.cosh(r/5)**(-1)), 0.05 * (np.cosh(r/5)**(-1)))),

        (1.58,
         (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/5)**(-1)), 0.05 * (np.cosh(r/5)**(-1))),
         (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/5)**(-1)), 0.05 * (np.cosh(r/5)**(-1)))),

        (1.59,
         (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/5)**(-1)), 0.05 * (np.cosh(r/5)**(-1))),
         (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/5)**(-1)), 0.05 * (np.cosh(r/5)**(-1)))),
         
        #Following functions need larger space to plot (higher R)
        # (1.595,
        #  (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/5)**(-1)), 0.05 * (np.cosh(r/5)**(-1))),
        #  (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/5)**(-1)), 0.05 * (np.cosh(r/5)**(-1)))),

        # (1.596,
        #  (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/5)**(-1)), 0.05 * (np.cosh(r/5)**(-1))),
        #  (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/5)**(-1)), 0.05 * (np.cosh(r/5)**(-1)))),

        # (1.597,
        #  (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/5)**(-1)), 0.05 * (np.cosh(r/5)**(-1))),
        #  (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/5)**(-1)), 0.05 * (np.cosh(r/5)**(-1)))),

        # (1.598,
        #  (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/10)**(-1)), 0.05 * (np.cosh(r/10)**(-1))),
        #  (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/10)**(-1)), 0.05 * (np.cosh(r/10)**(-1)))),

        # (1.599,
        #  (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/10)**(-1)), 0.05 * (np.cosh(r/10)**(-1))),
        #  (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/10)**(-1)), 0.05 * (np.cosh(r/10)**(-1)))),

        # (1.5995,
        #  (0.5 * np.ones(N+1) - 0.015 * (np.cosh(r/10)**(-1)), 0.025 * (np.cosh(r/10)**(-1))),
        #  (0.5 * np.ones(N+1) - 0.025 * (np.cosh(r/10)**(-1)), 0.025 * (np.cosh(r/10)**(-1)))),

        # (1.5997,
        #  (0.5 * np.ones(N+1) - 0.01 * (np.cosh(r/15)**(-1)), 0.01 * (np.cosh(r/15)**(-1))),
        #  (0.5 * np.ones(N+1) - 0.015 * (np.cosh(r/15)**(-1)), 0.015 * (np.cosh(r/15)**(-1)))),
         
        # (1.5998,
        #  (0.5 * np.ones(N+1) - 0.01 * (np.cosh(r/15)**(-1)), 0.01 * (np.cosh(r/15)**(-1))),
        #  (0.5 * np.ones(N+1) - 0.01 * (np.cosh(r/15)**(-1)), 0.01 * (np.cosh(r/15)**(-1)))),

        # (1.5999,
        #  (0.5 * np.ones(N+1) - 0.005 * (np.cosh(r/20)**(-1)), 0.005 * (np.cosh(r/20)**(-1))),
        #  (0.5 * np.ones(N+1) - 0.01 * (np.cosh(r/20)**(-1)), 0.01 * (np.cosh(r/20)**(-1)))),
    ]

    results_2d = []
    results_3d = []

    Q_values_2d = []
    E_values_2d = []
    Q_values_3d = []
    E_values_3d = []
    omega_list = []

    # Loop over all frequencies
    for omega, (f0_init_2d, f1_init_2d), (f0_init_3d, f1_init_3d) in omega_configs:
        print(f"\n{'='*60}")
        print(f"Processing ω = {omega}")
        print('='*60)

        # Solve in 2D
        r_2d, f0_2d, f1_2d = solve_2d_polar(
            omega=omega,
            f0_initial=f0_init_2d,
            f1_initial=f1_init_2d,
            R=R, N=N, mu=mu, lamb=lamb, lamb4=lamb4
        )
        Q_2d = compute_Q(r_2d, f1_2d, omega, dim=2)
        E_2d = compute_E(r_2d, f0_2d, f1_2d, omega, mu, lamb, lamb4, dim=2)
        results_2d.append((omega, r_2d, f0_2d, f1_2d, f0_init_2d, f1_init_2d, Q_2d, E_2d))
        Q_values_2d.append(Q_2d)
        E_values_2d.append(E_2d)

        # Solve in 3D
        r_3d, f0_3d, f1_3d = solve_3d_spherical(
            omega=omega,
            f0_initial=f0_init_3d,
            f1_initial=f1_init_3d,
            R=R, N=N, mu=mu, lamb=lamb, lamb4=lamb4
        )
        Q_3d = compute_Q(r_3d, f1_3d, omega, dim=3)
        E_3d = compute_E(r_3d, f0_3d, f1_3d, omega, mu, lamb, lamb4, dim=3)
        results_3d.append((omega, r_3d, f0_3d, f1_3d, f0_init_3d, f1_init_3d, Q_3d, E_3d))
        Q_values_3d.append(Q_3d)
        E_values_3d.append(E_3d)

        omega_list.append(omega)

    # Plot comparison of profiles (2D and 3D)
    print("\n=== Plotting profile solutions ===")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    colors = ['b', 'r', 'g', 'm', 'c', 'y', 'orange', 'purple', 'brown', 'pink']

    # 2D f0
    ax1 = axes[0, 0]
    for idx, (omega, r, f0, f1, _, _, Q, E) in enumerate(results_2d):
        ax1.plot(r, f0, color=colors[idx % len(colors)], linewidth=2, label=f'ω={omega:.2f}')
    ax1.set_title('2D Polar: f₀(r)')
    ax1.set_xlabel('r')
    ax1.set_ylabel('f₀(r)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # 2D f1
    ax2 = axes[0, 1]
    for idx, (omega, r, f0, f1, _, _, Q, E) in enumerate(results_2d):
        ax2.plot(r, f1, color=colors[idx % len(colors)], linewidth=2, label=f'ω={omega:.2f}')
    ax2.set_title('2D Polar: f₁(r)')
    ax2.set_xlabel('r')
    ax2.set_ylabel('f₁(r)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # 3D f0
    ax3 = axes[1, 0]
    for idx, (omega, r, f0, f1, _, _, Q, E) in enumerate(results_3d):
        ax3.plot(r, f0, color=colors[idx % len(colors)], linewidth=2, label=f'ω={omega:.2f}')
    ax3.set_title('3D Spherical: f₀(r)')
    ax3.set_xlabel('r')
    ax3.set_ylabel('f₀(r)')
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    # 3D f1
    ax4 = axes[1, 1]
    for idx, (omega, r, f0, f1, _, _, Q, E) in enumerate(results_3d):
        ax4.plot(r, f1, color=colors[idx % len(colors)], linewidth=2, label=f'ω={omega:.2f}')
    ax4.set_title('3D Spherical: f₁(r)')
    ax4.set_xlabel('r')
    ax4.set_ylabel('f₁(r)')
    ax4.grid(True, alpha=0.3)
    ax4.legend()

    plt.tight_layout()
    plt.savefig('solutions_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

    # Plot Q and E vs ω (log-log)
    print("\n=== Plotting Q and E vs ω ===")
    omega_arr = np.array(omega_list)
    Q_2d_arr = np.array(Q_values_2d)
    E_2d_arr = np.array(E_values_2d)
    Q_3d_arr = np.array(Q_values_3d)
    E_3d_arr = np.array(E_values_3d)

    sort_idx = np.argsort(omega_arr)
    omega_sorted = omega_arr[sort_idx]
    Q_2d_sorted = Q_2d_arr[sort_idx]
    E_2d_sorted = E_2d_arr[sort_idx]
    Q_3d_sorted = Q_3d_arr[sort_idx]
    E_3d_sorted = E_3d_arr[sort_idx]

    fig2, axes2 = plt.subplots(2, 2, figsize=(14, 10))
    # 2D Q vs ω
    axes2[0,0].plot(np.log(omega_sorted), np.log(Q_2d_sorted), 'b.-', linewidth=2, markersize=8)
    axes2[0,0].set_title('2D: log(Q) vs log(ω)')
    axes2[0,0].set_xlabel('log(ω)')
    axes2[0,0].set_ylabel('log(Q)')
    axes2[0,0].grid(True, alpha=0.3)

    # 2D E vs ω
    axes2[0,1].plot(np.log(omega_sorted), np.log(E_2d_sorted), 'r.-', linewidth=2, markersize=8)
    axes2[0,1].set_title('2D: log(E) vs log(ω)')
    axes2[0,1].set_xlabel('log(ω)')
    axes2[0,1].set_ylabel('log(E)')
    axes2[0,1].grid(True, alpha=0.3)

    # 3D Q vs ω
    axes2[1,0].plot(np.log(omega_sorted), np.log(Q_3d_sorted), 'b.-', linewidth=2, markersize=8)
    axes2[1,0].set_title('3D: log(Q) vs log(ω)')
    axes2[1,0].set_xlabel('log(ω)')
    axes2[1,0].set_ylabel('log(Q)')
    axes2[1,0].grid(True, alpha=0.3)

    # 3D E vs ω
    axes2[1,1].plot(np.log(omega_sorted), np.log(E_3d_sorted), 'r.-', linewidth=2, markersize=8)
    axes2[1,1].set_title('3D: log(E) vs log(ω)')
    axes2[1,1].set_xlabel('log(ω)')
    axes2[1,1].set_ylabel('log(E)')
    axes2[1,1].grid(True, alpha=0.3)

    plt.tight_layout()
    # plt.savefig('Q_E_vs_omega.png', dpi=300, bbox_inches='tight')
    plt.show()

    # Save data files
    print("\n=== Saving data ===")
    for omega, r, f0, f1, f0_init, f1_init, Q, E in results_2d:
        data_2d = np.column_stack((r, f0, f1, f0_init, f1_init))
        np.savetxt(f'2D_solutions_omega_{omega:.2f}_4fields.txt', data_2d,
                   header=f'r f0 f1 f0_initial f1_initial (ω={omega}) Q={Q:.6e} E={E:.6e}',
                   fmt='%.8f', comments='')
        print(f"Saved 2D data for ω={omega:.2f}, Q={Q:.6e}, E={E:.6e}")

    for omega, r, f0, f1, f0_init, f1_init, Q, E in results_3d:
        data_3d = np.column_stack((r, f0, f1, f0_init, f1_init))
        np.savetxt(f'3D_solutions_omega_{omega:.2f}_4fields.txt', data_3d,
                   header=f'r f0 f1 f0_initial f1_initial (ω={omega}) Q={Q:.6e} E={E:.6e}',
                   fmt='%.8f', comments='')
        print(f"Saved 3D data for ω={omega:.2f}, Q={Q:.6e}, E={E:.6e}")


if __name__ == "__main__":
    main()
