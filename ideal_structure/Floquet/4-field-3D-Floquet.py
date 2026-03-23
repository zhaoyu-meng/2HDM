#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
4-field-3D-Floquet.py

Floquet stability analysis for the 4-field oscillon solution in 3D spherical coordinates.
This script solves the background profiles f₀(r) and f₁(r) for the oscillon,
then constructs the linearized evolution operator and computes the monodromy matrix
via high‑precision time integration. Finally, it computes the Floquet multipliers
and plots their distribution to assess linear stability.

Corresponds to the 4‑field model (two complex scalars) described in the paper,
with parameters given in Appendix 1. The 3D analysis uses the spherical Laplacian
∇² = d²/dr² + (2/r) d/dr.

Author: Zhaoyu Meng
Date:   June 2025
Paper:  Long-Lived Oscillons as Closed Domain Walls in the Z₂‑Symmetric Two‑Higgs‑Doublet Model
        arXiv:2601.14676 [hep-ph]
License: MIT (see LICENSE file)
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
# 1. System Parameters and Background Solution (3D Spherical)
# =============================================================================
class OscillonSystem3D:
    """
    Defines the 4‑field oscillon system in 3D spherical coordinates.
    Contains parameters, radial grid, background fields f0(r) and f1(r),
    and the radial Laplacian operator.
    """

    def __init__(self):
        # Physical parameters (dimensionless, see Appendix 1)
        self.mu = 0.5          # mass parameter μ²
        self.lam = 1.78        # quartic coupling λ
        self.lam4 = -2.56      # coupling λ₄
        self.omega = 0.45      # oscillon frequency

        # Spatial grid (3D spherical coordinates)
        self.R_max = 32.0      # maximum radius [1/m_h]
        self.Nr = 129          # number of radial points
        self.r = np.linspace(0, self.R_max, self.Nr)
        self.dr = self.r[1] - self.r[0]   # radial step

        # Safe copy of r to avoid division by zero at origin
        self.r_safe = self.r.copy()
        self.r_safe[0] = 1e-10

        # Time period of the oscillon
        self.T = 2 * np.pi / self.omega

        # Number of fields (four real scalar fields)
        self.n_fields = 4
        self.total_dof = self.n_fields * self.Nr   # total spatial degrees of freedom

        # Solve the background equations
        self.f0, self.f1 = self._solve_background_fields_3d()
        self.alpha = 0.0       # constant phase (does not affect stability)

        # Construct the 3D spherical Laplacian matrix
        self.D2_3d = self._build_3d_laplacian()

    def _solve_background_fields_3d(self):
        """
        Solve the static coupled ODEs for f₀(r) and f₁(r) in 3D spherical coordinates.
        The equations are given by Eq. (4.5) of the paper, adapted for 3D.
        """
        print("Solving background fields in 3D spherical coordinates...")

        N = self.Nr
        dr = self.dr
        r = self.r
        r_safe = self.r_safe

        def equations_3d(y):
            f0 = y[:N]
            f1 = y[N:]
            res0 = np.zeros(N)
            res1 = np.zeros(N)

            # 3D spherical Laplacian: ∇²f = f'' + (2/r)f'
            for i in range(N):
                if i == 0:
                    # At r=0, use limit: ∇²f ≈ 3 * f''(0)
                    # with f''(0) approximated by 2*(f[1] - f[0])/dr²
                    laplacian_f0 = 6.0 * (f0[1] - f0[0]) / dr**2
                    laplacian_f1 = 6.0 * (f1[1] - f1[0]) / dr**2
                elif i == N - 1:
                    # Neumann boundary condition at r=R_max
                    laplacian_f0 = 2.0 * (f0[N-2] - f0[N-1]) / dr**2
                    laplacian_f1 = 2.0 * (f1[N-2] - f1[N-1]) / dr**2
                else:
                    # Finite differences for interior points
                    second_f0 = (f0[i-1] - 2*f0[i] + f0[i+1]) / dr**2
                    first_f0  = (f0[i+1] - f0[i-1]) / (2 * dr)
                    laplacian_f0 = second_f0 + (2 / r_safe[i]) * first_f0

                    second_f1 = (f1[i-1] - 2*f1[i] + f1[i+1]) / dr**2
                    first_f1  = (f1[i+1] - f1[i-1]) / (2 * dr)
                    laplacian_f1 = second_f1 + (2 / r_safe[i]) * first_f1

                # Nonlinear terms from the potential (same as in 2D)
                nonlinear_f0 = (self.mu * f0[i]
                                - 2 * self.lam * f0[i] * (2*f0[i]**2 + 2*f1[i]**2)
                                + 2 * self.lam4 * f0[i] * (-f0[i]**2 + f1[i]**2))
                nonlinear_f1 = (self.omega**2 * f1[i] + self.mu * f1[i]
                                - 2 * self.lam * f1[i] * (2*f0[i]**2 + 2*f1[i]**2)
                                - 2 * self.lam4 * f1[i] * (-f0[i]**2 + f1[i]**2))

                res0[i] = laplacian_f0 + nonlinear_f0
                res1[i] = laplacian_f1 + nonlinear_f1

            return np.concatenate([res0, res1])

        # Initial guess based on hyperbolic tangent profiles (tuned for 3D)
        f0_init = (1 + np.tanh((self.r - 8.5) / 2)) * 0.2
        f1_init = (1 - np.tanh((self.r - 8.5) / 2)) * 0.25
        y_init = np.concatenate([f0_init, f1_init])

        print("Solving nonlinear equations for f₀ and f₁ in 3D...")
        try:
            sol = fsolve(equations_3d, y_init, xtol=1e-11, maxfev=5_000_000)
            print("3D background solution found!")
        except Exception as e:
            print(f"Warning: fsolve failed with error: {e}")
            print("Using initial guess as solution")
            sol = y_init

        f0_sol = sol[:self.Nr]
        f1_sol = sol[self.Nr:]

        return f0_sol, f1_sol

    def _build_3d_laplacian(self):
        """
        Build the sparse matrix representing the 3D spherical Laplacian operator:
        ∇² = d²/dr² + (2/r) d/dr, with Neumann boundary conditions at r=0 and r=R_max.
        """
        print("Building 3D spherical Laplacian matrix...")
        N = self.Nr
        dr = self.dr
        rows, cols, data = [], [], []

        # Interior points (i=1 to N-2)
        for i in range(1, N - 1):
            r_i = self.r_safe[i]
            # Second derivative term
            rows.append(i); cols.append(i-1); data.append(1.0 / dr**2)
            rows.append(i); cols.append(i);   data.append(-2.0 / dr**2)
            rows.append(i); cols.append(i+1); data.append(1.0 / dr**2)

            # First derivative term (2/r) d/dr
            if r_i > 0:
                coeff = 1.0 / (r_i * dr)   # factor 2/(2dr) = 1/(r dr)
                rows.append(i); cols.append(i-1); data.append(-coeff)
                rows.append(i); cols.append(i+1); data.append( coeff)

        # Boundary at r=0 (Neumann)
        # ∇²f ≈ 6(f₁ - f₀)/dr²
        rows.append(0); cols.append(0); data.append(-6.0 / dr**2)
        rows.append(0); cols.append(1); data.append( 6.0 / dr**2)

        # Boundary at r=R_max (Neumann)
        rows.append(N-1); cols.append(N-2); data.append( 2.0 / dr**2)
        rows.append(N-1); cols.append(N-1); data.append(-2.0 / dr**2)

        D2 = sp.coo_matrix((data, (rows, cols)), shape=(N, N), dtype=np.complex128).tocsr()
        return D2

    def background_fields(self, t):
        """
        Return the background field configuration at time t.
        The fields are constructed from f₀(r) and f₁(r) using the ansatz
        given in Eq. (4.2) with α = 0.
        """
        cos_wt = np.cos(self.omega * t)
        sin_wt = np.sin(self.omega * t)
        phi = np.zeros((self.n_fields, self.Nr), dtype=np.complex128)
        phi[0, :] = self.f0 * np.cos(self.alpha) + self.f1 * cos_wt
        phi[2, :] = -self.f0 * np.sin(self.alpha) - self.f1 * sin_wt
        phi[1, :] = -self.f0 * np.cos(self.alpha) + self.f1 * cos_wt
        phi[3, :] = self.f0 * np.sin(self.alpha) - self.f1 * sin_wt
        return phi

    def compute_potential_jacobian(self, t):
        """
        Compute the Jacobian matrix J(t, r) = ∂²V/∂ψ_i∂ψ_j of the potential,
        evaluated on the background solution at time t. The result is a
        (n_fields, n_fields, Nr) array.
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

            # First field (ψ₁)
            J[0, 0, k] = mu - 2*lam*(3*p1**2 + p5**2 + s34) - 2*lam4*p2**2
            J[0, 1, k] = -4*lam*p1*p2 - 2*lam4*(c + p1*p2)
            J[0, 2, k] = -4*lam*p1*p5 - 2*lam4*p2*p6
            J[0, 3, k] = -4*lam*p1*p6 - 2*lam4*p2*p5

            # Second field (ψ₂)
            J[1, 0, k] = -4*lam*p1*p2 - 2*lam4*(c + p1*p2)
            J[1, 1, k] = mu - 2*lam*(3*p2**2 + p6**2 + s12) - 2*lam4*p1**2
            J[1, 2, k] = -4*lam*p2*p5 - 2*lam4*p1*p6
            J[1, 3, k] = -4*lam*p2*p6 - 2*lam4*p1*p5

            # Third field (ψ₅)
            J[2, 0, k] = -4*lam*p1*p5 - 2*lam4*p2*p6
            J[2, 1, k] = -4*lam*p2*p5 - 2*lam4*p1*p6
            J[2, 2, k] = mu - 2*lam*(3*p5**2 + p1**2 + s34) - 2*lam4*p6**2
            J[2, 3, k] = -4*lam*p5*p6 - 2*lam4*(c + p5*p6)

            # Fourth field (ψ₆)
            J[3, 0, k] = -4*lam*p1*p6 - 2*lam4*p2*p5
            J[3, 1, k] = -4*lam*p2*p6 - 2*lam4*p1*p5
            J[3, 2, k] = -4*lam*p5*p6 - 2*lam4*(c + p5*p6)
            J[3, 3, k] = mu - 2*lam*(3*p6**2 + p2**2 + s12) - 2*lam4*p5**2

        return J


# =============================================================================
# 2. Floquet Analyzer for 3D
# =============================================================================
class FloquetAnalyzer:
    """
    Performs Floquet stability analysis by constructing the linearized evolution
    operator and computing the monodromy matrix. The analysis includes a sponge
    layer to absorb outgoing perturbations.
    """

    def __init__(self, system):
        self.system = system
        self.Nr = system.Nr
        self.n_fields = system.n_fields
        self.total_dof = system.total_dof

    def build_linear_operator(self, t):
        """
        Build the full linear operator M(t) for the first‑order system
        d/dt (δ, δ̇)ᵀ = M(t) (δ, δ̇)ᵀ, where
        M = [[0, I], [D₂ + J(t), -Γ(r)]].
        The term -Γ(r) implements the sponge absorbing boundary condition.
        """
        # Laplacian part (block diagonal for each field)
        D_list = [self.system.D2_3d for _ in range(self.n_fields)]
        D_block = sp.block_diag(D_list, format='csr', dtype=np.complex128)

        # Jacobian part (potential curvature)
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

        # Sponge layer profile (quadratic ramp near the boundary)
        r = self.system.r
        R_max = self.system.R_max
        absorb_width = R_max * 0.1
        absorb_start = R_max * 0.9
        sponge_strength = 5.0
        sponge_profile = np.zeros_like(r)
        mask = r > absorb_start
        if np.any(mask):
            ratio = (r[mask] - absorb_start) / absorb_width
            sponge_profile[mask] = sponge_strength * ratio**2
        S_diag = sp.diags(sponge_profile, format='csr', dtype=np.complex128)
        S_block = sp.block_diag([S_diag] * self.n_fields, format='csr')

        # Assemble full 2N×2N matrix
        dim = self.total_dof
        Zero = sp.csr_matrix((dim, dim), dtype=np.complex128)
        Eye = sp.eye(dim, format='csr', dtype=np.complex128)
        Full_L = sp.bmat([[Zero, Eye],
                          [L_force, -S_block]], format='csr')
        return Full_L

    def compute_full_monodromy_matrix_high_precision(self, n_time_samples=50):
        """
        Compute the monodromy matrix M = U(T, 0) by integrating the linearized
        equations column by column with high precision. The integration uses
        the DOP853 method (high-order Runge‑Kutta) and directly evaluates
        the linear operator at each time step (no pre‑computed interpolation).
        """
        T = self.system.T
        dim = 2 * self.total_dof
        print(f"Computing full monodromy matrix with high precision (with Sponge Layer)...")
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
                                max_step=T/140, dense_output=False)
                M[:, col] = sol.y[:, -1]
            except Exception as e:
                print(f"  Warning: Integration failed for column {col}: {e}")
                print(f"  Using identity column as fallback")

        print("Monodromy matrix computation completed!")
        return M

    def compute_all_eigenvalues_high_precision(self):
        """
        Compute the monodromy matrix and then all its eigenvalues and eigenvectors.
        This is the most time‑consuming step; use with care.
        """
        print("\nComputing full monodromy matrix with high precision...")
        M = self.compute_full_monodromy_matrix_high_precision(n_time_samples=200000)

        print("\nComputing all eigenvalues and eigenvectors with high precision...")
        print("This may take a long time, please be patient...")
        M_dense = M.toarray() if sp.issparse(M) else M
        eigenvalues, eigenvectors = np.linalg.eig(M_dense)

        # Sort by magnitude (largest first)
        sorted_idx = np.argsort(-np.abs(eigenvalues))
        eigenvalues = eigenvalues[sorted_idx]
        eigenvectors = eigenvectors[:, sorted_idx]
        print(f"Computation completed! Found {len(eigenvalues)} eigenvalues and eigenvectors")

        return eigenvalues, eigenvectors, M

    def analyze_stability_high_precision(self):
        """
        Full stability analysis: compute eigenvalues, Floquet exponents,
        and classify modes as stable, marginal, or unstable.
        """
        eigenvalues, eigenvectors, M = self.compute_all_eigenvalues_high_precision()

        T = self.system.T
        exponents = np.log(eigenvalues) / T
        moduli = np.abs(eigenvalues)

        # Classification thresholds
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
            sorted_unstable_idx = unstable_indices[np.argsort(-unstable_growth_rates)]
            results['sorted_unstable_indices'] = sorted_unstable_idx
            results['most_unstable_index'] = sorted_unstable_idx[0]
        else:
            results['max_growth_rate'] = 0.0
            results['avg_growth_rate'] = 0.0
            results['max_modulus'] = np.max(moduli)
            results['sorted_unstable_indices'] = np.array([])
            results['most_unstable_index'] = None

        return results


# =============================================================================
# 3. Plotting Functions (Publication Quality)
# =============================================================================
def plot_most_unstable_mode_pattern(system, results):
    """
    Plot the spatial pattern of the most unstable mode (if any).
    This function visualizes the real part, imaginary part, magnitude, phase,
    energy distribution among fields, and a 2D polar projection of the dominant field.
    """
    if results['num_unstable'] == 0:
        print("\nNo unstable modes found, cannot plot unstable mode pattern.")
        return

    most_unstable_idx = results['most_unstable_index']
    eigenvector = results['eigenvectors'][:, most_unstable_idx]

    print(f"\nAnalyzing the most unstable mode (index {most_unstable_idx}):")
    print(f"  Growth rate: {np.real(results['exponents'][most_unstable_idx]):.6e}")
    print(f"  |ρ|: {np.abs(results['multipliers'][most_unstable_idx]):.6e}")

    # Extract position components (first half of eigenvector)
    total_dof = system.total_dof
    pos_components = eigenvector[:total_dof]

    # Split into the four fields
    field_names = ['ψ₁', 'ψ₂', 'ψ₅', 'ψ₆']
    field_patterns = []
    for i in range(system.n_fields):
        start = i * system.Nr
        end = (i + 1) * system.Nr
        field_patterns.append(pos_components[start:end])

    # Create comprehensive plot
    fig = plt.figure(figsize=(16, 12))

    # 1. Real parts
    ax1 = plt.subplot(3, 2, 1)
    for i, (name, pattern) in enumerate(zip(field_names, field_patterns)):
        ax1.plot(system.r, np.real(pattern), label=f'{name} (real)', linewidth=2, alpha=0.8)
    ax1.set_xlabel('Radial coordinate r', fontsize=12)
    ax1.set_ylabel('Amplitude (real part)', fontsize=12)
    ax1.set_title('Real Parts of Unstable Mode Patterns', fontsize=14)
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([0, system.R_max/2])

    # 2. Imaginary parts
    ax2 = plt.subplot(3, 2, 2)
    for i, (name, pattern) in enumerate(zip(field_names, field_patterns)):
        ax2.plot(system.r, np.imag(pattern), label=f'{name} (imag)', linewidth=2, alpha=0.8)
    ax2.set_xlabel('Radial coordinate r', fontsize=12)
    ax2.set_ylabel('Amplitude (imaginary part)', fontsize=12)
    ax2.set_title('Imaginary Parts of Unstable Mode Patterns', fontsize=14)
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([0, system.R_max/2])

    # 3. Magnitudes
    ax3 = plt.subplot(3, 2, 3)
    for i, (name, pattern) in enumerate(zip(field_names, field_patterns)):
        ax3.plot(system.r, np.abs(pattern), label=f'{name}', linewidth=2, alpha=0.8)
    ax3.set_xlabel('Radial coordinate r', fontsize=12)
    ax3.set_ylabel('Magnitude', fontsize=12)
    ax3.set_title('Magnitude of Unstable Mode Patterns', fontsize=14)
    ax3.legend(loc='upper right', fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim([0, system.R_max/2])

    # 4. Phases (unwrapped)
    ax4 = plt.subplot(3, 2, 4)
    for i, (name, pattern) in enumerate(zip(field_names, field_patterns)):
        phase = np.angle(pattern)
        phase_unwrapped = np.unwrap(phase)
        ax4.plot(system.r, phase_unwrapped, label=f'{name}', linewidth=2, alpha=0.8)
    ax4.set_xlabel('Radial coordinate r', fontsize=12)
    ax4.set_ylabel('Phase (radians)', fontsize=12)
    ax4.set_title('Phase of Unstable Mode Patterns', fontsize=14)
    ax4.legend(loc='upper right', fontsize=10)
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim([0, system.R_max/2])

    # 5. Energy distribution among fields
    ax5 = plt.subplot(3, 2, 5)
    field_energies = [np.sum(np.abs(p)**2) for p in field_patterns]
    total_energy = np.sum(field_energies)
    energy_ratios = [e / total_energy for e in field_energies]
    bars = ax5.bar(field_names, energy_ratios,
                   color=['blue', 'red', 'green', 'purple'], alpha=0.7)
    ax5.set_xlabel('Field', fontsize=12)
    ax5.set_ylabel('Energy Ratio', fontsize=12)
    ax5.set_title('Energy Distribution Among Fields', fontsize=14)
    ax5.grid(True, alpha=0.3, axis='y')
    for bar, ratio in zip(bars, energy_ratios):
        height = bar.get_height()
        ax5.text(bar.get_x() + bar.get_width()/2., height,
                 f'{ratio:.3f}', ha='center', va='bottom', fontsize=10)

    # 6. 2D polar projection of the dominant field
    ax6 = plt.subplot(3, 2, 6)
    dominant_idx = np.argmax(energy_ratios)
    dominant_pattern = field_patterns[dominant_idx]
    theta = np.linspace(0, 2*np.pi, 100)
    r = system.r
    # Use only first half of radial points for clarity
    R, Theta = np.meshgrid(r[:len(r)//2], theta)
    pattern_mag = np.abs(dominant_pattern[:len(r)//2])
    pattern_2d = np.outer(np.ones_like(theta), pattern_mag)
    im = ax6.contourf(Theta, R, pattern_2d, levels=50, cmap='viridis')
    ax6.set_aspect('equal')
    ax6.set_title(f'2D Pattern of Dominant Field ({field_names[dominant_idx]})', fontsize=14)
    ax6.set_xlabel('θ', fontsize=12)
    ax6.set_ylabel('r', fontsize=12)
    plt.colorbar(im, ax=ax6, shrink=0.8, label='Amplitude')

    plt.suptitle(f'Most Unstable Mode Analysis (Mode {most_unstable_idx})', fontsize=16, y=0.98)
    plt.tight_layout()
    plt.savefig('most_unstable_mode_pattern_4field.png', dpi=300, bbox_inches='tight')
    plt.show()

    # Print detailed info
    print(f"\nMost Unstable Mode Analysis:")
    print(f"  Mode index: {most_unstable_idx}")
    print(f"  Growth rate: {np.real(results['exponents'][most_unstable_idx]):.6e}")
    print(f"  Frequency: {np.imag(results['exponents'][most_unstable_idx]):.6e}")
    print(f"  |ρ|: {np.abs(results['multipliers'][most_unstable_idx]):.6e}")
    print(f"  arg(ρ): {np.angle(results['multipliers'][most_unstable_idx]):.6f} rad")
    print(f"\nEnergy distribution among fields:")
    for name, ratio in zip(field_names, energy_ratios):
        print(f"  {name}: {ratio:.3f} ({ratio*100:.1f}%)")
    print(f"\nDominant field: {field_names[dominant_idx]} ({energy_ratios[dominant_idx]*100:.1f}% of total energy)")


def plot_background_fields_3d(system):
    """Plot the background profiles f₀(r) and f₁(r) in 3D."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.plot(system.r, system.f0, 'b-', linewidth=2)
    ax1.set_xlabel('Radial coordinate r', fontsize=12)
    ax1.set_ylabel('Field amplitude f₀(r)', fontsize=12)
    ax1.set_title('Background field f₀(r) in 3D', fontsize=14)
    ax1.grid(True, alpha=0.3)

    ax2.plot(system.r, system.f1, 'r-', linewidth=2)
    ax2.set_xlabel('Radial coordinate r', fontsize=12)
    ax2.set_ylabel('Field amplitude f₁(r)', fontsize=12)
    ax2.set_title('Background field f₁(r) in 3D', fontsize=14)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('background_fields_3d_4field.png', dpi=150, bbox_inches='tight')
    plt.show()


def create_paper_format_plots_3d(system, results):
    """
    Create publication‑quality plots for 3D spherical coordinates:
    - Left: Floquet multipliers on the unit circle.
    - Right: Histogram of |ρ|.
    """
    eigenvalues = results['multipliers']
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Unit circle
    theta = np.linspace(0, 2*np.pi, 200)
    ax1.plot(np.cos(theta), np.sin(theta), 'k-', linewidth=1.0, alpha=0.7, label='Unit circle')

    moduli = np.abs(eigenvalues)
    stable_mask = moduli <= 0.99999
    one_mask = (moduli > 0.99999) & (moduli < 1.00001)
    unstable_mask = moduli >= 1.00001

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
    ax1.set_title('Floquet Multipliers (3D Spherical)', fontsize=16)
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.axis('equal')
    ax1.set_xlim([-1.5, 1.5])
    ax1.set_ylim([-1.5, 1.5])

    # Histogram
    bins = np.linspace(0, max(2, np.max(moduli)*1.1), 30)
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

    stats_text = f'Total eigenvalues: {len(eigenvalues)}\n'
    stats_text += f'Stable (|ρ| < 1): {len(stable_eig)}\n'
    stats_text += f'Marginal (|ρ| = 1): {len(one_eig)}\n'
    stats_text += f'Unstable (|ρ| > 1): {len(unstable_eig)}'
    ax2.text(0.65, 0.95, stats_text, transform=ax2.transAxes,
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))

    plt.tight_layout()
    plt.savefig('floquet_analysis_3d_simple_4field.png', dpi=300)
    plt.show()

    print("\nAnalysis completed:")
    print(f"  Total eigenvalues analyzed: {len(eigenvalues)}")
    print(f"  Stable eigenvalues (|ρ| < 1): {len(stable_eig)}")
    print(f"  Marginal eigenvalues (|ρ| = 1): {len(one_eig)}")
    print(f"  Unstable eigenvalues (|ρ| > 1): {len(unstable_eig)}")
    if len(unstable_eig) > 0:
        print(f"  Maximum modulus: {np.max(moduli):.6f}")


# =============================================================================
# 4. Main Function
# =============================================================================
def main_3d():
    print("=" * 60)
    print("Floquet Stability Analysis in 3D Spherical Coordinates")
    print("Parameters: μ=0.5, λ=1.78, λ₄=-2.56, ω=0.45")
    print("=" * 60)

    print("\nInitializing system in 3D spherical coordinates...")
    system = OscillonSystem3D()

    print("\nPlotting background fields...")
    plot_background_fields_3d(system)

    print("\nCreating Floquet analyzer...")
    analyzer = FloquetAnalyzer(system)

    print("\nCalculating eigenvalues...")
    results = analyzer.analyze_stability_high_precision()

    print("\nGenerating publication‑quality plots...")
    create_paper_format_plots_3d(system, results)

    plot_most_unstable_mode_pattern(system, results)


if __name__ == "__main__":
    main_3d()
