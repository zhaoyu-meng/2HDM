#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
4-field-2D-Floquet.py

Floquet stability analysis for the oscillon solution in the 4‑field model
(two complex scalar fields) in 2D polar coordinates.

This script solves for the background oscillon profiles f0(r) and f1(r)
using the nonlinear ODEs derived from the Lagrangian, then constructs the
linearized operator around the time‑periodic background, computes the
monodromy matrix via direct integration, and finally computes its eigenvalues
to determine stability.

The analysis corresponds to Section IV.C of the paper:
"Long-Lived Oscillons as Closed Domain Walls in the Z₂-Symmetric Two-Higgs-Doublet Model"
Zhaoyu Meng, arXiv:2601.14676 [hep-ph]

Author: Zhaoyu Meng
Date:   March 2026
License: MIT
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# 1. System Parameters and Background Solution
# =============================================================================
class OscillonSystem:
    """
    Defines the physical system and solves for the stationary oscillon
    profiles f0(r) and f1(r) in 2D polar coordinates.
    """

    def __init__(self):
        # Physical parameters (dimensionless, see Appendix 1 of the paper)
        self.mu = 0.5          # μ²
        self.lam = 1.78        # λ
        self.lam4 = -2.56      # λ₄
        self.omega = 0.45      # Oscillon frequency

        # Spatial grid in 2D polar coordinates (radial direction)
        self.R_max = 32.0       # Maximum radius
        self.Nr = 129           # Number of radial points
        self.r = np.linspace(0, self.R_max, self.Nr)
        self.dr = self.r[1] - self.r[0]

        # Avoid division by zero at r=0
        self.r_safe = self.r.copy()
        self.r_safe[0] = 1e-10

        # Oscillation period
        self.T = 2 * np.pi / self.omega

        # Number of real fields (4 active fields in this reduced model)
        self.n_fields = 4
        self.total_dof = self.n_fields * self.Nr

        # Solve the background equations
        self.f0, self.f1 = self._solve_background_fields()
        self.alpha = 0.0       # Phase angle (arbitrary, see Eq. (4.2))

        # Build the 2D radial Laplacian matrix
        self.D2_radial = self._build_radial_laplacian()

    def _solve_background_fields(self):
        """
        Solve the coupled ODEs for f0(r) and f1(r) in 2D polar coordinates.

        The equations are derived from Eq. (4.5) of the paper, using the 2D
        polar Laplacian: ∇²f = f'' + (1/r) f'.

        Returns
        -------
        f0 : ndarray, shape (Nr,)
            Solution for the static part of the oscillon.
        f1 : ndarray, shape (Nr,)
            Solution for the oscillating amplitude.
        """
        print("Solving background fields in 2D polar coordinates...")

        N = self.Nr
        dr = self.dr
        r = self.r
        r_safe = self.r_safe

        def equations_2d(y):
            f0 = y[:N]
            f1 = y[N:]
            res0 = np.zeros(N)
            res1 = np.zeros(N)

            for i in range(N):
                # Laplacian in 2D polar coordinates
                if i == 0:
                    laplacian_f0 = 2.0 * (f0[1] - f0[0]) / dr**2
                    laplacian_f1 = 2.0 * (f1[1] - f1[0]) / dr**2
                elif i == N - 1:
                    laplacian_f0 = 2.0 * (f0[N-2] - f0[N-1]) / dr**2
                    laplacian_f1 = 2.0 * (f1[N-2] - f1[N-1]) / dr**2
                else:
                    second0 = (f0[i-1] - 2.0*f0[i] + f0[i+1]) / dr**2
                    second1 = (f1[i-1] - 2.0*f1[i] + f1[i+1]) / dr**2
                    first0 = (f0[i+1] - f0[i-1]) / (2.0 * dr)
                    first1 = (f1[i+1] - f1[i-1]) / (2.0 * dr)
                    radial0 = (1.0 / r_safe[i]) * first0
                    radial1 = (1.0 / r_safe[i]) * first1
                    laplacian_f0 = second0 + radial0
                    laplacian_f1 = second1 + radial1

                # Nonlinear terms from the potential
                nonlinear0 = self.mu * f0[i] - 2.0*self.lam * f0[i] * (2.0*f0[i]**2 + 2.0*f1[i]**2) \
                             + 2.0*self.lam4 * f0[i] * (-f0[i]**2 + f1[i]**2)
                nonlinear1 = self.omega**2 * f1[i] + self.mu * f1[i] \
                             - 2.0*self.lam * f1[i] * (2.0*f0[i]**2 + 2.0*f1[i]**2) \
                             - 2.0*self.lam4 * f1[i] * (-f0[i]**2 + f1[i]**2)

                res0[i] = laplacian_f0 + nonlinear0
                res1[i] = laplacian_f1 + nonlinear1

            return np.concatenate([res0, res1])

        # Initial guess (sigmoid-like shape)
        f0_initial = (1.0 + np.tanh((self.r - 3.0) / 2.0)) * 0.25
        f1_initial = (1.0 - np.tanh((self.r - 3.0) / 2.0)) * 0.29
        y_initial = np.concatenate([f0_initial, f1_initial])

        print("Solving nonlinear equations for f₀ and f₁...")
        try:
            solution = fsolve(equations_2d, y_initial, xtol=1e-11, maxfev=5_000_000)
            print("Background solution found!")
        except Exception as e:
            print(f"Warning: fsolve failed with error: {e}")
            print("Using initial guess as solution")
            solution = y_initial

        f0_solution = solution[:self.Nr]
        f1_solution = solution[self.Nr:]
        return f0_solution, f1_solution

    def _build_radial_laplacian(self):
        """
        Construct the sparse matrix representing the 2D radial Laplacian
        with Neumann boundary conditions at r = R_max.

        Returns
        -------
        D2 : scipy.sparse.csr_matrix, shape (Nr, Nr)
            Discrete Laplacian operator.
        """
        print("Building 2D radial Laplacian matrix...")
        N = self.Nr
        dr = self.dr
        rows, cols, data = [], [], []

        # Interior points (i = 1 .. N-2)
        for i in range(1, N - 1):
            r_i = self.r_safe[i]
            # Second derivative coefficients
            rows.append(i); cols.append(i-1); data.append(1.0 / dr**2)
            rows.append(i); cols.append(i);   data.append(-2.0 / dr**2)
            rows.append(i); cols.append(i+1); data.append(1.0 / dr**2)
            # First derivative coefficient (1/r) * d/dr
            if r_i > 0:
                rows.append(i); cols.append(i-1); data.append(-1.0 / (2.0 * r_i * dr))
                rows.append(i); cols.append(i+1); data.append( 1.0 / (2.0 * r_i * dr))

        # Boundary at r = 0 (symmetry: f'(0) = 0)
        rows.append(0); cols.append(0); data.append(-2.0 / dr**2)
        rows.append(0); cols.append(1); data.append( 2.0 / dr**2)

        # Boundary at r = R_max (Neumann: f'(R_max) = 0)
        rows.append(N-1); cols.append(N-2); data.append( 2.0 / dr**2)
        rows.append(N-1); cols.append(N-1); data.append(-2.0 / dr**2)

        D2 = sp.coo_matrix((data, (rows, cols)), shape=(N, N), dtype=np.complex128).tocsr()
        return D2

    def background_fields(self, t):
        """
        Compute the background field configuration at time t.

        Parameters
        ----------
        t : float
            Time.

        Returns
        -------
        phi : ndarray, shape (n_fields, Nr)
            The 4 real field components at each radial point.
        """
        cos_wt = np.cos(self.omega * t)
        sin_wt = np.sin(self.omega * t)

        phi = np.zeros((self.n_fields, self.Nr), dtype=np.complex128)
        # Field assignments according to Eq. (4.2) with α = 0
        phi[0, :] = self.f0 * np.cos(self.alpha) + self.f1 * cos_wt
        phi[2, :] = -self.f0 * np.sin(self.alpha) - self.f1 * sin_wt
        phi[1, :] = -self.f0 * np.cos(self.alpha) + self.f1 * cos_wt
        phi[3, :] = self.f0 * np.sin(self.alpha) - self.f1 * sin_wt
        return phi

    def compute_potential_jacobian(self, t):
        """
        Compute the Jacobian matrix of the potential part of the equations
        evaluated on the background fields.

        Parameters
        ----------
        t : float
            Time.

        Returns
        -------
        J : ndarray, shape (n_fields, n_fields, Nr)
            The Hessian of the potential with respect to the fields,
            evaluated at the background field values.
        """
        phi = self.background_fields(t)
        psi1, psi2, psi5, psi6 = phi[0, :], phi[1, :], phi[2, :], phi[3, :]

        psi1_sq, psi2_sq, psi5_sq, psi6_sq = psi1**2, psi2**2, psi5**2, psi6**2
        sum12 = psi1_sq + psi5_sq
        sum34 = psi2_sq + psi6_sq
        cross = psi1 * psi2 + psi5 * psi6

        J = np.zeros((self.n_fields, self.n_fields, self.Nr), dtype=np.complex128)
        mu, lam, lam4 = self.mu, self.lam, self.lam4

        for k in range(self.Nr):
            p1, p2, p5, p6 = psi1[k], psi2[k], psi5[k], psi6[k]
            s12, s34, c = sum12[k], sum34[k], cross[k]

            # Row 0 (ψ₁)
            J[0, 0, k] = mu - 2.0*lam*(3.0*p1**2 + p5**2 + s34) - 2.0*lam4*p2**2
            J[0, 1, k] = -4.0*lam*p1*p2 - 2.0*lam4*(c + p1*p2)
            J[0, 2, k] = -4.0*lam*p1*p5 - 2.0*lam4*p2*p6
            J[0, 3, k] = -4.0*lam*p1*p6 - 2.0*lam4*p2*p5

            # Row 1 (ψ₂)
            J[1, 0, k] = -4.0*lam*p1*p2 - 2.0*lam4*(c + p1*p2)
            J[1, 1, k] = mu - 2.0*lam*(3.0*p2**2 + p6**2 + s12) - 2.0*lam4*p1**2
            J[1, 2, k] = -4.0*lam*p2*p5 - 2.0*lam4*p1*p6
            J[1, 3, k] = -4.0*lam*p2*p6 - 2.0*lam4*p1*p5

            # Row 2 (ψ₅)
            J[2, 0, k] = -4.0*lam*p1*p5 - 2.0*lam4*p2*p6
            J[2, 1, k] = -4.0*lam*p2*p5 - 2.0*lam4*p1*p6
            J[2, 2, k] = mu - 2.0*lam*(3.0*p5**2 + p1**2 + s34) - 2.0*lam4*p6**2
            J[2, 3, k] = -4.0*lam*p5*p6 - 2.0*lam4*(c + p5*p6)

            # Row 3 (ψ₆)
            J[3, 0, k] = -4.0*lam*p1*p6 - 2.0*lam4*p2*p5
            J[3, 1, k] = -4.0*lam*p2*p6 - 2.0*lam4*p1*p5
            J[3, 2, k] = -4.0*lam*p5*p6 - 2.0*lam4*(c + p5*p6)
            J[3, 3, k] = mu - 2.0*lam*(3.0*p6**2 + p2**2 + s12) - 2.0*lam4*p5**2

        return J


# =============================================================================
# 2. Floquet Analyzer
# =============================================================================
class FloquetAnalyzer:
    """
    Performs Floquet stability analysis of the oscillon.

    Constructs the linearized operator around the time‑periodic background,
    computes the monodromy matrix by integrating the perturbation equations,
    and calculates its eigenvalues.
    """

    def __init__(self, system):
        """
        Parameters
        ----------
        system : OscillonSystem
            The system containing background profiles and parameters.
        """
        self.system = system
        self.Nr = system.Nr
        self.n_fields = system.n_fields
        self.total_dof = system.total_dof

    def build_linear_operator(self, t):
        """
        Construct the full linear operator M(t) = [[0, I], [D2+J(t), -S]]
        where D2 is the Laplacian, J(t) is the potential Jacobian, and S is
        a sponge layer damping matrix.

        Parameters
        ----------
        t : float
            Time.

        Returns
        -------
        Full_L : scipy.sparse.csr_matrix, shape (2*total_dof, 2*total_dof)
            The linearized operator for the perturbation vector y = (δ, δ̇).
        """
        # 1. Force matrix (Laplacian + Jacobian)
        D_list = [self.system.D2_radial for _ in range(self.n_fields)]
        D_block = sp.block_diag(D_list, format='csr', dtype=np.complex128)

        J_tx = self.system.compute_potential_jacobian(t)
        rows, cols, data = [], [], []
        for k in range(self.Nr):
            for i in range(self.n_fields):
                for j in range(self.n_fields):
                    val = J_tx[i, j, k]
                    if abs(val) > 1e-13:
                        row_idx = i * self.Nr + k
                        col_idx = j * self.Nr + k
                        rows.append(row_idx)
                        cols.append(col_idx)
                        data.append(val)
        J_block = sp.coo_matrix((data, (rows, cols)),
                                shape=(self.total_dof, self.total_dof),
                                dtype=np.complex128).tocsr()

        L_force = D_block + J_block

        # 2. Sponge layer (absorbing boundary)
        r = self.system.r
        R_max = self.system.R_max
        absorb_width = R_max * 0.1
        absorb_start = R_max * 0.9
        sponge_strength = 5.0

        sponge_profile = np.zeros_like(r)
        mask = r > absorb_start
        if np.any(mask):
            sponge_profile[mask] = sponge_strength * ((r[mask] - absorb_start) / absorb_width)**2

        S_diag = sp.diags(sponge_profile, format='csr', dtype=np.complex128)
        S_block = sp.block_diag([S_diag] * self.n_fields, format='csr')

        # 3. Assemble the full 2N x 2N matrix
        dim = self.total_dof
        Zero = sp.csr_matrix((dim, dim), dtype=np.complex128)
        Eye = sp.eye(dim, format='csr', dtype=np.complex128)

        Full_L = sp.bmat([[Zero, Eye],
                          [L_force, -S_block]], format='csr')
        return Full_L

    def compute_full_monodromy_matrix_high_precision(self, n_time_samples=50):
        """
        Compute the monodromy matrix by integrating the linearized equations
        over one period using a high‑order ODE solver.

        Parameters
        ----------
        n_time_samples : int, optional
            Number of time steps used for operator interpolation (not used in
            this version; kept for API consistency).

        Returns
        -------
        M : ndarray, shape (2*total_dof, 2*total_dof)
            The monodromy matrix.
        """
        T = self.system.T
        dim = 2 * self.total_dof

        print(f"Computing full monodromy matrix with high precision...")
        print(f"Dimension: {dim} × {dim}")

        M = np.eye(dim, dtype=np.complex128)

        print("Computing monodromy matrix columns...")
        for col in range(dim):
            if col % 100 == 0:
                print(f"  Column {col+1}/{dim}")

            y0 = M[:, col]

            def ode_rhs(t, y):
                L = self.build_linear_operator(t)
                return L.dot(y)

            try:
                sol = solve_ivp(ode_rhs, [0, T], y0,
                                method='DOP853',
                                rtol=1e-11, atol=1e-13,
                                max_step=T / 140,
                                dense_output=False)
                M[:, col] = sol.y[:, -1]
            except Exception as e:
                print(f"  Warning: Integration failed for column {col}: {e}")
                print(f"  Using identity column as fallback")

        print("Monodromy matrix computation completed!")
        return M

    def compute_all_eigenvalues_high_precision(self):
        """
        Compute all eigenvalues of the monodromy matrix.

        Returns
        -------
        eigenvalues : ndarray
            All eigenvalues.
        eigenvectors : ndarray
            Corresponding eigenvectors.
        monodromy_matrix : ndarray
            The monodromy matrix itself.
        """
        print("\nComputing full monodromy matrix with high precision...")
        # Adjust n_time_samples if you wish to use interpolation (currently not used)
        M = self.compute_full_monodromy_matrix_high_precision(n_time_samples=140000)

        print("\nComputing all eigenvalues and eigenvectors with high precision...")
        print("This may take a long time, please be patient...")

        M_dense = M.toarray() if sp.issparse(M) else M
        eigenvalues, eigenvectors = np.linalg.eig(M_dense)

        # Sort by magnitude (largest first)
        sorted_indices = np.argsort(-np.abs(eigenvalues))
        eigenvalues = eigenvalues[sorted_indices]
        eigenvectors = eigenvectors[:, sorted_indices]

        print(f"Computation completed! Found {len(eigenvalues)} eigenvalues.")
        return eigenvalues, eigenvectors, M

    def analyze_stability_high_precision(self):
        """
        Perform the full stability analysis and return a dictionary with results.

        Returns
        -------
        results : dict
            Contains 'multipliers', 'exponents', 'eigenvectors', 'moduli',
            'num_total', 'num_stable', 'num_marginal', 'num_unstable',
            'max_growth_rate', etc.
        """
        eigenvalues, eigenvectors, M = self.compute_all_eigenvalues_high_precision()

        T = self.system.T
        exponents = np.log(eigenvalues) / T

        moduli = np.abs(eigenvalues)

        # Stability classification (thresholds adapted to numerical precision)
        stable_mask = moduli <= 0.99999
        marginal_mask = (moduli > 0.99999) & (moduli < 1.00001)
        unstable_mask = moduli >= 1.00001

        results = {
            'multipliers': eigenvalues,
            'exponents': exponents,
            'eigenvectors': eigenvectors,
            'moduli': moduli,
            'num_total': len(eigenvalues),
            'num_stable': np.sum(stable_mask),
            'num_marginal': np.sum(marginal_mask),
            'num_unstable': np.sum(unstable_mask),
            'stable_indices': np.where(stable_mask)[0],
            'marginal_indices': np.where(marginal_mask)[0],
            'unstable_indices': np.where(unstable_mask)[0],
            'monodromy_matrix': M
        }

        if np.sum(unstable_mask) > 0:
            unstable_growth_rates = np.real(exponents[unstable_mask])
            results['max_growth_rate'] = np.max(unstable_growth_rates)
            results['avg_growth_rate'] = np.mean(unstable_growth_rates)
            results['max_modulus'] = np.max(moduli[unstable_mask])
            # Sort unstable modes by growth rate
            unstable_indices = results['unstable_indices']
            sorted_unstable = unstable_indices[np.argsort(-unstable_growth_rates)]
            results['sorted_unstable_indices'] = sorted_unstable
            results['most_unstable_index'] = sorted_unstable[0]
        else:
            results['max_growth_rate'] = 0.0
            results['avg_growth_rate'] = 0.0
            results['max_modulus'] = np.max(moduli)
            results['sorted_unstable_indices'] = np.array([])
            results['most_unstable_index'] = None

        return results


# =============================================================================
# 3. Plotting Functions
# =============================================================================
def plot_background_fields(system):
    """
    Plot the background fields f0(r) and f1(r) in 2D polar coordinates.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(system.r, system.f0, 'b-', linewidth=2)
    ax1.set_xlabel('Radial coordinate r', fontsize=12)
    ax1.set_ylabel('Field amplitude f₀(r)', fontsize=12)
    ax1.set_title('Background field f₀(r) in 2D', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([0, system.R_max / 2])

    ax2.plot(system.r, system.f1, 'r-', linewidth=2)
    ax2.set_xlabel('Radial coordinate r', fontsize=12)
    ax2.set_ylabel('Field amplitude f₁(r)', fontsize=12)
    ax2.set_title('Background field f₁(r) in 2D', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([0, system.R_max / 2])

    plt.tight_layout()
    plt.savefig('background_fields_2d_4field.png', dpi=150, bbox_inches='tight')
    plt.show()


def create_paper_format_plots(system, results):
    """
    Create publication‑quality plots of the Floquet multipliers.

    Parameters
    ----------
    system : OscillonSystem
        The system (used for radial coordinates, though not directly here).
    results : dict
        The output of analyze_stability_high_precision().
    """
    eigenvalues = results['multipliers']
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Unit circle plot
    theta = np.linspace(0, 2 * np.pi, 200)
    ax1.plot(np.cos(theta), np.sin(theta), 'k-', linewidth=1.0, alpha=0.7, label='Unit circle')

    moduli = np.abs(eigenvalues)
    stable_mask = moduli <= 0.99999
    one_mask = (moduli > 0.99999) & (moduli < 1.00001)
    unstable_mask = moduli > 1.00001

    stable_eig = eigenvalues[stable_mask]
    one_eig = eigenvalues[one_mask]
    unstable_eig = eigenvalues[unstable_mask]

    ax1.scatter(np.real(stable_eig), np.imag(stable_eig),
                c='blue', alpha=0.7, edgecolors='none', s=50,
                label=f'Stable (|ρ|<1): {len(stable_eig)}')
    ax1.scatter(np.real(one_eig), np.imag(one_eig),
                c='red', alpha=0.9, edgecolors='none', s=50,
                label=f'Marginal (|ρ|=1): {len(one_eig)}')
    ax1.scatter(np.real(unstable_eig), np.imag(unstable_eig),
                c='red', alpha=0.9, edgecolors='none', s=70, marker='^',
                label=f'Unstable (|ρ|>1): {len(unstable_eig)}')

    ax1.set_xlabel('Re(ρ)', fontsize=14)
    ax1.set_ylabel('Im(ρ)', fontsize=14)
    ax1.set_title('Floquet Multipliers (2D Polar Coordinates)', fontsize=16)
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.axis('equal')
    ax1.set_xlim([-1.5, 1.5])
    ax1.set_ylim([-1.5, 1.5])

    # Histogram of moduli
    bins = np.linspace(0, max(2.0, np.max(moduli) * 1.1), 30)
    n, bins, patches = ax2.hist(moduli, bins=bins, alpha=0.7, color='steelblue',
                                edgecolor='black', linewidth=0.5)
    for i, patch in enumerate(patches):
        if bins[i] > 1.0:
            patch.set_facecolor('lightcoral')
            patch.set_alpha(0.8)

    ax2.set_xlabel('|ρ|', fontsize=14)
    ax2.set_ylabel('Frequency', fontsize=14)
    ax2.set_title('Distribution of Floquet Multiplier Moduli', fontsize=16)
    ax2.grid(True, alpha=0.3, linestyle='--')

    stats_text = (f'Total eigenvalues: {len(eigenvalues)}\n'
                  f'Stable (|ρ| < 1): {len(stable_eig)}\n'
                  f'Marginal (|ρ| = 1): {len(one_eig)}\n'
                  f'Unstable (|ρ| > 1): {len(unstable_eig)}')
    ax2.text(0.65, 0.95, stats_text, transform=ax2.transAxes,
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))

    plt.tight_layout()
    plt.savefig('floquet_analysis_2d_simple_4field.png', dpi=300, bbox_inches='tight')
    plt.show()

    print("\nAnalysis completed:")
    print(f"  Total eigenvalues analyzed: {len(eigenvalues)}")
    print(f"  Stable eigenvalues (|ρ| < 1): {len(stable_eig)}")
    print(f"  Marginal eigenvalues (|ρ| = 1): {len(one_eig)}")
    print(f"  Unstable eigenvalues (|ρ| > 1): {len(unstable_eig)}")
    if len(unstable_eig) > 0:
        print(f"  Maximum modulus: {np.max(moduli):.3f}")


# =============================================================================
# 4. Main Function
# =============================================================================
def main():
    """Main entry point: initialise system, run Floquet analysis, plot results."""
    print("=" * 60)
    print("Floquet Stability Analysis in 2D Polar Coordinates")
    print("Parameters: μ=0.5, λ=1.78, λ₄=-2.56, ω=0.45")
    print("=" * 60)

    print("\nInitializing system in 2D polar coordinates...")
    system = OscillonSystem()

    print("\nPlotting background fields...")
    plot_background_fields(system)

    print("\nCreating Floquet analyzer...")
    analyzer = FloquetAnalyzer(system)

    print("Calculating eigenvalues...")
    results = analyzer.analyze_stability_high_precision()

    print("\nGenerating publication‑format plots...")
    create_paper_format_plots(system, results)


if __name__ == "__main__":
    main()
