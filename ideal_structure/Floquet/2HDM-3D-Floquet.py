#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2HDM-3D-Floquet.py

Floquet stability analysis for the full 2HDM oscillon solution in 3D spherical coordinates.
This script solves the background profiles f₀(r) and f₁(r) for the oscillon in the
Z₂-symmetric two-Higgs-doublet model (2HDM) with 8 real scalar fields. It then
constructs the linearized evolution operator (including a sponge absorbing layer)
and computes the monodromy matrix via high‑precision time integration. The Floquet
multipliers are computed and plotted to assess linear stability.

The parameters correspond to the radiation‑free choice discussed in the paper,
with λ₁ = λ₂ = λ₃/2, λ₄ = -|λ₅|, and all masses set to 200 GeV (see Appendix 2).

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
class EightFieldOscillonSystem3D:
    """
    Defines the full 2HDM oscillon system in 3D spherical coordinates.
    Contains parameters, radial grid, background fields f0(r) and f1(r),
    and the radial Laplacian operator.
    """

    def __init__(self):
        # Physical parameters (dimensionless, see Appendix 2)
        self.mu1 = 0.5        # μ₁²
        self.mu2 = 0.5        # μ₂² (equal in the radiation‑free case)

        # Coupling constants satisfying λ₁ = λ₂ = λ₃/2 and λ₄ + |λ₅| = 0
        self.lam3 = 3.56      # λ₃
        self.lam1 = self.lam3 / 2.0
        self.lam2 = self.lam3 / 2.0
        self.lam4 = -2.56     # λ₄
        self.lam5_abs = 2.56  # |λ₅|

        self.omega = 0.45     # oscillon frequency (in units of m_h)

        # Spatial grid (3D spherical coordinates)
        self.R_max = 32.0     # maximum radius [1/m_h]
        self.Nr = 129         # number of radial points
        self.r = np.linspace(0, self.R_max, self.Nr)
        self.dr = self.r[1] - self.r[0]

        # Avoid division by zero at r=0
        self.r_safe = self.r.copy()
        self.r_safe[0] = 1e-10

        # Time period
        self.T = 2 * np.pi / self.omega

        # Number of fields (8 real scalar fields)
        self.n_fields = 8
        self.total_dof = self.n_fields * self.Nr

        # Solve background equations
        self.f0, self.f1 = self._solve_background_fields_3d()

        # Phase angles (see Eq. (5.3) of the paper)
        self.alpha1 = 0.5     # α (first doublet)
        self.alpha2 = 0.2     # α (second doublet)
        self.beta = 0.4       # β (mixing angle between doublets)
        self.theta_alpha = 0.7  # θ_β (additional phase)

        # Construct 3D spherical Laplacian matrix
        self.D2_3d = self._build_3d_spherical_laplacian()

        # Verify background solution
        self._verify_background()

    def _solve_background_fields_3d(self):
        """
        Solve the static coupled ODEs for f₀(r) and f₁(r) in 3D spherical coordinates.
        The equations are derived from Eq. (5.5) of the paper, adapted for 3D.
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
                    # Interior points: central differences
                    second_f0 = (f0[i-1] - 2*f0[i] + f0[i+1]) / dr**2
                    first_f0  = (f0[i+1] - f0[i-1]) / (2 * dr)
                    laplacian_f0 = second_f0 + (2 / r_safe[i]) * first_f0

                    second_f1 = (f1[i-1] - 2*f1[i] + f1[i+1]) / dr**2
                    first_f1  = (f1[i+1] - f1[i-1]) / (2 * dr)
                    laplacian_f1 = second_f1 + (2 / r_safe[i]) * first_f1

                # Nonlinear terms (Eq. (5.5))
                nonlinear_f0 = (self.mu1 * f0[i]
                                - self.lam3 * f0[i] * (2*f0[i]**2 + 2*f1[i]**2)
                                + 2 * self.lam4 * f0[i] * (-f0[i]**2 + f1[i]**2))
                nonlinear_f1 = (self.omega**2 * f1[i] + self.mu1 * f1[i]
                                - self.lam3 * f1[i] * (2*f0[i]**2 + 2*f1[i]**2)
                                - 2 * self.lam4 * f1[i] * (-f0[i]**2 + f1[i]**2))

                res0[i] = laplacian_f0 + nonlinear_f0
                res1[i] = laplacian_f1 + nonlinear_f1

            return np.concatenate([res0, res1])

        # Initial guess based on hyperbolic tangent profiles (tuned for 3D)
        f0_init = (1 + np.tanh((self.r - 4.5) / 2)) * 0.2
        f1_init = (1 - np.tanh((self.r - 4.5) / 2)) * 0.25
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

    def _build_3d_spherical_laplacian(self):
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

    def _verify_background(self):
        """Print basic verification of the background solution."""
        print("\nBackground solution verification:")
        print(f"  f₀ range: [{np.min(self.f0):.6f}, {np.max(self.f0):.6f}]")
        print(f"  f₁ range: [{np.min(self.f1):.6f}, {np.max(self.f1):.6f}]")
        print(f"  Center: f₀(0)={self.f0[0]:.6f}, f₁(0)={self.f1[0]:.6f}")
        print(f"  Boundary: f₀(R)={self.f0[-1]:.6f}, f₁(R)={self.f1[-1]:.6f}")

    def background_fields(self, t):
        """
        Return the background field configuration at time t.
        Constructed from f₀(r) and f₁(r) using the ansatz in Eq. (5.3) of the paper.
        """
        phi = np.zeros((self.n_fields, self.Nr), dtype=np.complex128)

        # First doublet (components 0-3)
        phi[0, :] = (+self.f0 * np.cos(self.alpha1) + self.f1 * np.cos(self.omega * t)) * np.cos(self.beta)
        phi[1, :] = (-self.f0 * np.cos(self.alpha1) + self.f1 * np.cos(self.omega * t)) * np.cos(self.beta)
        phi[2, :] = (+self.f0 * np.cos(self.alpha2) - self.f1 * np.cos(self.omega * t + self.theta_alpha)) * np.sin(self.beta)
        phi[3, :] = (-self.f0 * np.cos(self.alpha2) - self.f1 * np.cos(self.omega * t + self.theta_alpha)) * np.sin(self.beta)

        # Second doublet (components 4-7)
        phi[4, :] = (-self.f0 * np.sin(self.alpha1) - self.f1 * np.sin(self.omega * t)) * np.cos(self.beta)
        phi[5, :] = (+self.f0 * np.sin(self.alpha1) - self.f1 * np.sin(self.omega * t)) * np.cos(self.beta)
        phi[6, :] = (+self.f0 * np.sin(self.alpha2) + self.f1 * np.sin(self.omega * t + self.theta_alpha)) * np.sin(self.beta)
        phi[7, :] = (-self.f0 * np.sin(self.alpha2) + self.f1 * np.sin(self.omega * t + self.theta_alpha)) * np.sin(self.beta)

        return phi

    def compute_potential_jacobian(self, t):
        """
        Compute the Jacobian matrix J(t, r) = ∂²V/∂ψ_i∂ψ_j of the potential,
        evaluated on the background solution at time t.
        Returns an array of shape (n_fields, n_fields, Nr).
        """
        phi = self.background_fields(t)
        psi = phi   # psi[0] = ψ₁, ..., psi[7] = ψ₈

        # Parameters
        mu1, mu2 = self.mu1, self.mu2
        lam1, lam2, lam3 = self.lam1, self.lam2, self.lam3
        lam4 = self.lam4

        # Since λ₄ + |λ₅| = 0, we have λ₄ - |λ₅| = 2λ₄
        lam4_minus_lam5 = 2 * lam4

        J = np.zeros((self.n_fields, self.n_fields, self.Nr), dtype=np.complex128)

        for k in range(self.Nr):
            p = psi[:, k]   # field values at this radial point

            # Useful combinations
            sum1 = p[0]**2 + p[2]**2 + p[4]**2 + p[6]**2   # Ψ₁†Ψ₁
            sum2 = p[1]**2 + p[3]**2 + p[5]**2 + p[7]**2   # Ψ₂†Ψ₂
            re_cross = p[0]*p[1] + p[4]*p[5] + p[2]*p[3] + p[6]*p[7]  # Re(Ψ₁†Ψ₂)

            # Jacobian elements (∂²V/∂ψ_i∂ψ_j)
            # Component 0 (ψ₁)
            J[0, 0, k] = mu1 - 2*lam1*(3*p[0]**2 + p[2]**2 + p[4]**2 + p[6]**2 + sum1) - lam3*sum2 - lam4_minus_lam5*p[1]**2
            J[0, 1, k] = -4*lam1*p[0]*p[1] - lam4_minus_lam5*(re_cross + p[0]*p[1])
            J[0, 2, k] = -4*lam1*p[0]*p[2] - lam4_minus_lam5*p[1]*p[3]
            J[0, 3, k] = -4*lam1*p[0]*p[3] - lam4_minus_lam5*p[1]*p[2]
            J[0, 4, k] = -4*lam1*p[0]*p[4] - lam4_minus_lam5*p[1]*p[5]
            J[0, 5, k] = -4*lam1*p[0]*p[5] - lam4_minus_lam5*p[1]*p[4]
            J[0, 6, k] = -4*lam1*p[0]*p[6] - lam4_minus_lam5*p[1]*p[7]
            J[0, 7, k] = -4*lam1*p[0]*p[7] - lam4_minus_lam5*p[1]*p[6]

            # Component 1 (ψ₂)
            J[1, 0, k] = -4*lam2*p[1]*p[0] - lam4_minus_lam5*(re_cross + p[0]*p[1])
            J[1, 1, k] = mu2 - 2*lam2*(3*p[1]**2 + p[3]**2 + p[5]**2 + p[7]**2 + sum2) - lam3*sum1 - lam4_minus_lam5*p[0]**2
            J[1, 2, k] = -4*lam2*p[1]*p[2] - lam4_minus_lam5*p[0]*p[3]
            J[1, 3, k] = -4*lam2*p[1]*p[3] - lam4_minus_lam5*p[0]*p[2]
            J[1, 4, k] = -4*lam2*p[1]*p[4] - lam4_minus_lam5*p[0]*p[5]
            J[1, 5, k] = -4*lam2*p[1]*p[5] - lam4_minus_lam5*p[0]*p[4]
            J[1, 6, k] = -4*lam2*p[1]*p[6] - lam4_minus_lam5*p[0]*p[7]
            J[1, 7, k] = -4*lam2*p[1]*p[7] - lam4_minus_lam5*p[0]*p[6]

            # Component 2 (ψ₃)
            J[2, 0, k] = -4*lam1*p[2]*p[0] - lam4_minus_lam5*p[3]*p[1]
            J[2, 1, k] = -4*lam1*p[2]*p[1] - lam4_minus_lam5*p[3]*p[0]
            J[2, 2, k] = mu1 - 2*lam1*(p[0]**2 + 3*p[2]**2 + p[4]**2 + p[6]**2 + sum1) - lam3*sum2 - lam4_minus_lam5*p[3]**2
            J[2, 3, k] = -4*lam1*p[2]*p[3] - lam4_minus_lam5*(re_cross + p[2]*p[3])
            J[2, 4, k] = -4*lam1*p[2]*p[4] - lam4_minus_lam5*p[3]*p[5]
            J[2, 5, k] = -4*lam1*p[2]*p[5] - lam4_minus_lam5*p[3]*p[4]
            J[2, 6, k] = -4*lam1*p[2]*p[6] - lam4_minus_lam5*p[3]*p[7]
            J[2, 7, k] = -4*lam1*p[2]*p[7] - lam4_minus_lam5*p[3]*p[6]

            # Component 3 (ψ₄)
            J[3, 0, k] = -4*lam2*p[3]*p[0] - lam4_minus_lam5*p[2]*p[1]
            J[3, 1, k] = -4*lam2*p[3]*p[1] - lam4_minus_lam5*p[2]*p[0]
            J[3, 2, k] = -4*lam2*p[3]*p[2] - lam4_minus_lam5*(re_cross + p[2]*p[3])
            J[3, 3, k] = mu2 - 2*lam2*(p[1]**2 + 3*p[3]**2 + p[5]**2 + p[7]**2 + sum2) - lam3*sum1 - lam4_minus_lam5*p[2]**2
            J[3, 4, k] = -4*lam2*p[3]*p[4] - lam4_minus_lam5*p[2]*p[5]
            J[3, 5, k] = -4*lam2*p[3]*p[5] - lam4_minus_lam5*p[2]*p[4]
            J[3, 6, k] = -4*lam2*p[3]*p[6] - lam4_minus_lam5*p[2]*p[7]
            J[3, 7, k] = -4*lam2*p[3]*p[7] - lam4_minus_lam5*p[2]*p[6]

            # Component 4 (ψ₅)
            J[4, 0, k] = -4*lam1*p[4]*p[0] - lam4_minus_lam5*p[5]*p[1]
            J[4, 1, k] = -4*lam1*p[4]*p[1] - lam4_minus_lam5*p[5]*p[0]
            J[4, 2, k] = -4*lam1*p[4]*p[2] - lam4_minus_lam5*p[5]*p[3]
            J[4, 3, k] = -4*lam1*p[4]*p[3] - lam4_minus_lam5*p[5]*p[2]
            J[4, 4, k] = mu1 - 2*lam1*(p[0]**2 + p[2]**2 + 3*p[4]**2 + p[6]**2 + sum1) - lam3*sum2 - lam4_minus_lam5*p[5]**2
            J[4, 5, k] = -4*lam1*p[4]*p[5] - lam4_minus_lam5*(re_cross + p[4]*p[5])
            J[4, 6, k] = -4*lam1*p[4]*p[6] - lam4_minus_lam5*p[5]*p[7]
            J[4, 7, k] = -4*lam1*p[4]*p[7] - lam4_minus_lam5*p[5]*p[6]

            # Component 5 (ψ₆)
            J[5, 0, k] = -4*lam2*p[5]*p[0] - lam4_minus_lam5*p[4]*p[1]
            J[5, 1, k] = -4*lam2*p[5]*p[1] - lam4_minus_lam5*p[4]*p[0]
            J[5, 2, k] = -4*lam2*p[5]*p[2] - lam4_minus_lam5*p[4]*p[3]
            J[5, 3, k] = -4*lam2*p[5]*p[3] - lam4_minus_lam5*p[4]*p[2]
            J[5, 4, k] = -4*lam2*p[5]*p[4] - lam4_minus_lam5*(re_cross + p[4]*p[5])
            J[5, 5, k] = mu2 - 2*lam2*(p[1]**2 + p[3]**2 + 3*p[5]**2 + p[7]**2 + sum2) - lam3*sum1 - lam4_minus_lam5*p[4]**2
            J[5, 6, k] = -4*lam2*p[5]*p[6] - lam4_minus_lam5*p[4]*p[7]
            J[5, 7, k] = -4*lam2*p[5]*p[7] - lam4_minus_lam5*p[4]*p[6]

            # Component 6 (ψ₇)
            J[6, 0, k] = -4*lam1*p[6]*p[0] - lam4_minus_lam5*p[7]*p[1]
            J[6, 1, k] = -4*lam1*p[6]*p[1] - lam4_minus_lam5*p[7]*p[0]
            J[6, 2, k] = -4*lam1*p[6]*p[2] - lam4_minus_lam5*p[7]*p[3]
            J[6, 3, k] = -4*lam1*p[6]*p[3] - lam4_minus_lam5*p[7]*p[2]
            J[6, 4, k] = -4*lam1*p[6]*p[4] - lam4_minus_lam5*p[7]*p[5]
            J[6, 5, k] = -4*lam1*p[6]*p[5] - lam4_minus_lam5*p[7]*p[4]
            J[6, 6, k] = mu1 - 2*lam1*(p[0]**2 + p[2]**2 + p[4]**2 + 3*p[6]**2 + sum1) - lam3*sum2 - lam4_minus_lam5*p[7]**2
            J[6, 7, k] = -4*lam1*p[6]*p[7] - lam4_minus_lam5*(re_cross + p[6]*p[7])

            # Component 7 (ψ₈)
            J[7, 0, k] = -4*lam2*p[7]*p[0] - lam4_minus_lam5*p[6]*p[1]
            J[7, 1, k] = -4*lam2*p[7]*p[1] - lam4_minus_lam5*p[6]*p[0]
            J[7, 2, k] = -4*lam2*p[7]*p[2] - lam4_minus_lam5*p[6]*p[3]
            J[7, 3, k] = -4*lam2*p[7]*p[3] - lam4_minus_lam5*p[6]*p[2]
            J[7, 4, k] = -4*lam2*p[7]*p[4] - lam4_minus_lam5*p[6]*p[5]
            J[7, 5, k] = -4*lam2*p[7]*p[5] - lam4_minus_lam5*p[6]*p[4]
            J[7, 6, k] = -4*lam2*p[7]*p[6] - lam4_minus_lam5*(re_cross + p[6]*p[7])
            J[7, 7, k] = mu2 - 2*lam2*(p[1]**2 + p[3]**2 + p[5]**2 + 3*p[7]**2 + sum2) - lam3*sum1 - lam4_minus_lam5*p[6]**2

        return J


# =============================================================================
# 2. Floquet Analyzer for 8-Field System (3D)
# =============================================================================
class EightFieldFloquetAnalyzer3D:
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
        equations column by column with high precision. Uses a fixed number of
        pre‑computed time samples and linear interpolation between them.
        """
        T = self.system.T
        dim = 2 * self.total_dof
        print(f"Computing full monodromy matrix with high precision (with Sponge Layer)...")
        print(f"Dimension: {dim} × {dim}")
        print(f"Using {n_time_samples} time samples")

        # Precompute linear operators at different times
        times = np.linspace(0, T, n_time_samples, endpoint=False)
        L_operators = []
        print("Precomputing linear operators at time samples...")
        for i, t in enumerate(times):
            if i % 1000 == 0:
                print(f"  Time sample {i+1}/{n_time_samples}: t = {t:.4f}")
            L_operators.append(self.build_linear_operator(t))

        # Initialize monodromy matrix as identity
        M = np.eye(dim, dtype=np.complex128)

        # Compute all columns
        print("Computing monodromy matrix columns...")
        for col in range(dim):
            if col % 100 == 0:
                print(f"  Column {col+1}/{dim}")
            y0 = M[:, col]

            def ode_rhs(t, y):
                # Find nearest time sample with linear interpolation
                idx = int(t / T * n_time_samples) % n_time_samples
                idx_next = (idx + 1) % n_time_samples
                t1 = times[idx]
                t2 = times[idx_next] if idx_next != 0 else T
                if t2 > t1:
                    alpha = (t - t1) / (t2 - t1)
                else:
                    alpha = (t - t1) / (T - t1 + t2)
                L1 = L_operators[idx]
                L2 = L_operators[idx_next] if idx_next < n_time_samples else L_operators[0]
                L = (1.0 - alpha) * L1 + alpha * L2
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
        M = self.compute_full_monodromy_matrix_high_precision(n_time_samples=40000)

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

        # Classification thresholds (adjusted for numerical precision)
        stable_mask = moduli <= 0.9999
        marginal_mask = (moduli > 0.9999) & (moduli < 1.0001)
        unstable_mask = moduli >= 1.0001

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
# 3. Plot Functions (Publication Quality for 3D)
# =============================================================================
def plot_most_unstable_mode_pattern_3d(system, results):
    """
    Visualize the spatial pattern of the most unstable Floquet mode (3D version).
    Creates six subplots: real/imag/magnitude/phase profiles, energy distribution,
    and a 2D polar contour of the dominant field component.
    """
    if results['num_unstable'] == 0:
        print("\nNo unstable modes found, cannot plot unstable mode pattern.")
        return

    # Get the most unstable mode
    most_unstable_idx = results['most_unstable_index']
    eigenvector = results['eigenvectors'][:, most_unstable_idx]

    print(f"\nAnalyzing the most unstable mode (index {most_unstable_idx}):")
    print(f"  Growth rate: {np.real(results['exponents'][most_unstable_idx]):.6e}")
    print(f"  |ρ|: {np.abs(results['multipliers'][most_unstable_idx]):.6e}")

    # Extract position components (first half of eigenvector)
    total_dof = system.total_dof
    pos_components = eigenvector[:total_dof]

    # Split into the 8 fields
    field_names = ['ψ1', 'ψ2', 'ψ3', 'ψ4', 'ψ5', 'ψ6', 'ψ7', 'ψ8']
    field_patterns = []
    for i in range(system.n_fields):
        start_idx = i * system.Nr
        end_idx = (i + 1) * system.Nr
        field_patterns.append(pos_components[start_idx:end_idx])

    # Create comprehensive plot
    fig = plt.figure(figsize=(16, 12))

    # Plot 1: Real parts
    ax1 = plt.subplot(3, 2, 1)
    for i, (name, pattern) in enumerate(zip(field_names, field_patterns)):
        ax1.plot(system.r, np.real(pattern), label=f'{name} (real)', linewidth=2, alpha=0.8)
    ax1.set_xlabel('Radial coordinate r', fontsize=12)
    ax1.set_ylabel('Amplitude (real part)', fontsize=12)
    ax1.set_title('Real Parts of Unstable Mode Patterns (3D)', fontsize=14)
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([0, system.R_max/2])

    # Plot 2: Imaginary parts
    ax2 = plt.subplot(3, 2, 2)
    for i, (name, pattern) in enumerate(zip(field_names, field_patterns)):
        ax2.plot(system.r, np.imag(pattern), label=f'{name} (imag)', linewidth=2, alpha=0.8)
    ax2.set_xlabel('Radial coordinate r', fontsize=12)
    ax2.set_ylabel('Amplitude (imaginary part)', fontsize=12)
    ax2.set_title('Imaginary Parts of Unstable Mode Patterns (3D)', fontsize=14)
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([0, system.R_max/2])

    # Plot 3: Magnitudes
    ax3 = plt.subplot(3, 2, 3)
    for i, (name, pattern) in enumerate(zip(field_names, field_patterns)):
        ax3.plot(system.r, np.abs(pattern), label=f'{name}', linewidth=2, alpha=0.8)
    ax3.set_xlabel('Radial coordinate r', fontsize=12)
    ax3.set_ylabel('Magnitude', fontsize=12)
    ax3.set_title('Magnitude of Unstable Mode Patterns (3D)', fontsize=14)
    ax3.legend(loc='upper right', fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim([0, system.R_max/2])

    # Plot 4: Phases
    ax4 = plt.subplot(3, 2, 4)
    for i, (name, pattern) in enumerate(zip(field_names, field_patterns)):
        phase = np.angle(pattern)
        phase_unwrapped = np.unwrap(phase)
        ax4.plot(system.r, phase_unwrapped, label=f'{name}', linewidth=2, alpha=0.8)
    ax4.set_xlabel('Radial coordinate r', fontsize=12)
    ax4.set_ylabel('Phase (radians)', fontsize=12)
    ax4.set_title('Phase of Unstable Mode Patterns (3D)', fontsize=14)
    ax4.legend(loc='upper right', fontsize=10)
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim([0, system.R_max/2])

    # Plot 5: Energy distribution among fields
    ax5 = plt.subplot(3, 2, 5)
    field_energies = [np.sum(np.abs(pattern)**2) for pattern in field_patterns]
    total_energy = np.sum(field_energies)
    energy_ratios = [e / total_energy for e in field_energies]

    bars = ax5.bar(field_names, energy_ratios, color=['blue', 'red', 'green', 'purple']*2, alpha=0.7)
    ax5.set_xlabel('Field', fontsize=12)
    ax5.set_ylabel('Energy Ratio', fontsize=12)
    ax5.set_title('Energy Distribution Among Fields', fontsize=14)
    ax5.grid(True, alpha=0.3, axis='y')
    for bar, ratio in zip(bars, energy_ratios):
        height = bar.get_height()
        ax5.text(bar.get_x() + bar.get_width()/2., height,
                 f'{ratio:.3f}', ha='center', va='bottom', fontsize=10)

    # Plot 6: 2D polar contour of the dominant field
    ax6 = plt.subplot(3, 2, 6)
    dominant_field_idx = np.argmax(energy_ratios)
    dominant_pattern = field_patterns[dominant_field_idx]

    theta = np.linspace(0, 2*np.pi, 100)
    r = system.r
    R, Theta = np.meshgrid(r[:len(r)//2], theta)   # only first half of radial points
    pattern_magnitude = np.abs(dominant_pattern[:len(r)//2])
    pattern_2d = np.outer(np.ones_like(theta), pattern_magnitude)

    im = ax6.contourf(Theta, R, pattern_2d, levels=50, cmap='viridis')
    ax6.set_aspect('equal')
    ax6.set_title(f'2D Pattern of Dominant Field ({field_names[dominant_field_idx]})', fontsize=14)
    ax6.set_xlabel('θ', fontsize=12)
    ax6.set_ylabel('r', fontsize=12)
    plt.colorbar(im, ax=ax6, shrink=0.8, label='Amplitude')

    plt.suptitle(f'Most Unstable Mode Analysis (Mode {most_unstable_idx})', fontsize=16, y=0.98)
    plt.tight_layout()
    plt.savefig('most_unstable_mode_pattern_2hdm_3d.png', dpi=300, bbox_inches='tight')
    plt.show()

    # Print summary
    print(f"\nMost Unstable Mode Analysis:")
    print(f"  Mode index: {most_unstable_idx}")
    print(f"  Growth rate: {np.real(results['exponents'][most_unstable_idx]):.6e}")
    print(f"  Frequency: {np.imag(results['exponents'][most_unstable_idx]):.6e}")
    print(f"  |ρ|: {np.abs(results['multipliers'][most_unstable_idx]):.6e}")
    print(f"  arg(ρ): {np.angle(results['multipliers'][most_unstable_idx]):.6f} rad")
    print(f"\nEnergy distribution among fields:")
    for name, ratio in zip(field_names, energy_ratios):
        print(f"  {name}: {ratio:.3f} ({ratio*100:.1f}%)")
    print(f"\nDominant field: {field_names[dominant_field_idx]} ({energy_ratios[dominant_field_idx]*100:.1f}% of total energy)")


def plot_background_fields_3d(system):
    """Plot the background profiles f₀(r) and f₁(r) for the 2HDM in 3D."""
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
    plt.savefig('background_fields_3d_2hdm.png', dpi=150, bbox_inches='tight')
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
    stable_mask = moduli <= 0.9999
    one_mask = (moduli > 0.9999) & (moduli < 1.0001)
    unstable_mask = moduli >= 1.0001

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
    plt.savefig('floquet_analysis_3d_simple_2hdm.png', dpi=300)
    plt.show()

    print("\nAnalysis completed:")
    print(f"  Total eigenvalues analyzed: {len(eigenvalues)}")
    print(f"  Stable eigenvalues (|ρ| < 1): {len(stable_eig)}")
    print(f"  Marginal eigenvalues (|ρ| = 1): {len(one_eig)}")
    print(f"  Unstable eigenvalues (|ρ| > 1): {len(unstable_eig)}")
    if len(unstable_eig) > 0:
        print(f"  Maximum modulus: {np.max(moduli):.3f}")


# =============================================================================
# 4. Main Program
# =============================================================================
def main():
    print("=" * 70)
    print("EIGHT-FIELD OSCILLON SYSTEM - FLOQUET STABILITY ANALYSIS (3D)")
    print("=" * 70)

    print("\n1. Initializing 3D system and solving background...")
    system = EightFieldOscillonSystem3D()

    print("\n2. Creating Floquet analyzer for 3D system...")
    analyzer = EightFieldFloquetAnalyzer3D(system)

    print("\n3. Computing stability...")
    results = analyzer.analyze_stability_high_precision()

    print("\n4. Generating publication‑quality plots...")
    create_paper_format_plots_3d(system, results)

    print("\n5. Plotting most unstable mode pattern...")
    plot_most_unstable_mode_pattern_3d(system, results)

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
