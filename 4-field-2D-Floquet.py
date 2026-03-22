# import numpy as np
# import scipy.sparse as sp
# import scipy.sparse.linalg as spla
# import matplotlib.pyplot as plt
# from scipy.integrate import solve_ivp
# from scipy.optimize import fsolve
# import warnings
# warnings.filterwarnings('ignore')

# # ================================
# # 1. System Parameters and Background Solution
# # ================================
# class OscillonSystem:
#     def __init__(self):
#         # Physical parameters
#         self.mu = 0.5
#         self.lam = 1.78
#         self.lam4 = -2.56
#         self.omega = 0.45
        
#         # Spatial parameters in 2D polar coordinates
#         self.R_max = 32.0  # Maximum radius
#         self.Nr = 321      # Number of radial points
#         self.r = np.linspace(0, self.R_max, self.Nr)
#         self.dr = self.r[1] - self.r[0]
        
#         # Avoid division by zero at r=0
#         self.r_safe = self.r.copy()
#         self.r_safe[0] = 1e-10
        
#         # Time period
#         self.T = 2 * np.pi / self.omega
        
#         # Number of fields
#         self.n_fields = 4
        
#         # Total degrees of freedom
#         self.total_dof = self.n_fields * self.Nr
        
#         # Solve for background fields
#         self.f0, self.f1 = self._solve_background_fields()
#         self.alpha = 0.0
        
#         # Construct 2D radial Laplacian matrix
#         self.D2_radial = self._build_radial_laplacian()
    
#     def _solve_background_fields(self):
#         """Solve for background fields f0 and f1 in 2D polar coordinates"""
#         print("Solving background fields in 2D polar coordinates...")
        
#         N = self.Nr
#         dr = self.dr
#         r = self.r
#         r_safe = self.r_safe
        
#         def equations_2d(y):
#             f0 = y[:N]
#             f1 = y[N:]
#             res0 = np.zeros(N)
#             res1 = np.zeros(N)
            
#             # 2D radial Laplacian: ∇²f = f'' + (1/r)f'
#             for i in range(N):
#                 # Boundary conditions at r=0
#                 if i == 0:
#                     # At r=0: f'(0) = 0, so ∇²f = 2 * (f1 - f0) / dr²
#                     laplacian_f0 = 2 * (f0[1] - f0[0]) / dr**2
#                     laplacian_f1 = 2 * (f1[1] - f1[0]) / dr**2
#                 # Boundary at r=R_max
#                 elif i == N-1:
#                     # Neumann boundary: f'(R_max) = 0
#                     laplacian_f0 = 2 * (f0[N-2] - f0[N-1]) / dr**2
#                     laplacian_f1 = 2 * (f1[N-2] - f1[N-1]) / dr**2
#                 else:
#                     # Interior points
#                     # f'' ≈ (f_{i-1} - 2f_i + f_{i+1}) / dr²
#                     second_deriv_f0 = (f0[i-1] - 2*f0[i] + f0[i+1]) / dr**2
#                     second_deriv_f1 = (f1[i-1] - 2*f1[i] + f1[i+1]) / dr**2
                    
#                     # (1/r)f' ≈ (1/r_i) * (f_{i+1} - f_{i-1})/(2dr)
#                     first_deriv_f0 = (f0[i+1] - f0[i-1]) / (2 * dr)
#                     first_deriv_f1 = (f1[i+1] - f1[i-1]) / (2 * dr)
                    
#                     radial_term_f0 = (1 / r_safe[i]) * first_deriv_f0
#                     radial_term_f1 = (1 / r_safe[i]) * first_deriv_f1
                    
#                     laplacian_f0 = second_deriv_f0 + radial_term_f0
#                     laplacian_f1 = second_deriv_f1 + radial_term_f1
                
#                 # Nonlinear terms
#                 nonlinear_f0 = self.mu * f0[i] - 2*self.lam * f0[i] * (2 * f0[i]**2 + 2 * f1[i]**2) + 2 * self.lam4 * f0[i] * (-f0[i]**2 + f1[i]**2)
#                 nonlinear_f1 = self.omega**2 * f1[i] + self.mu * f1[i] - 2*self.lam * f1[i] * (2 * f0[i]**2 + 2 * f1[i]**2) - 2 * self.lam4 * f1[i] * (-f0[i]**2 + f1[i]**2)
                
#                 res0[i] = laplacian_f0 + nonlinear_f0
#                 res1[i] = laplacian_f1 + nonlinear_f1
            
#             return np.concatenate([res0, res1])
        
#         # Initial guess - Gaussian profiles
        
#         f0_initial = (1+np.tanh((self.r-3)/2))*0.25
#         f1_initial = (1-np.tanh((self.r-3)/2))*0.29
#         y_initial = np.concatenate([f0_initial, f1_initial])
        
#         print("Solving nonlinear equations for f₀ and f₁...")
#         try:
#             solution = fsolve(equations_2d, y_initial, xtol=1e-10, maxfev=50000)
#             print("Background solution found!")
#         except Exception as e:
#             print(f"Warning: fsolve failed with error: {e}")
#             print("Using initial guess as solution")
#             solution = y_initial
        
#         # Extract solutions
#         f0_solution = solution[:self.Nr]
#         f1_solution = solution[self.Nr:]
        

        
#         return f0_solution, f1_solution
    
#     def _build_radial_laplacian(self):
#         """Build 2D radial Laplacian matrix: ∇² = d²/dr² + (1/r) d/dr"""
#         print("Building 2D radial Laplacian matrix...")
        
#         N = self.Nr
#         dr = self.dr
        
#         # Initialize sparse matrix components
#         rows, cols, data = [], [], []
        
#         # Interior points (i=1 to N-2)
#         for i in range(1, N-1):
#             r_i = self.r_safe[i]
            
#             # Second derivative: f'' ≈ (f_{i-1} - 2f_i + f_{i+1})/dr²
#             # Coefficient for f_{i-1}
#             rows.append(i); cols.append(i-1); data.append(1.0 / dr**2)
#             # Coefficient for f_i
#             rows.append(i); cols.append(i); data.append(-2.0 / dr**2)
#             # Coefficient for f_{i+1}
#             rows.append(i); cols.append(i+1); data.append(1.0 / dr**2)
            
#             # First derivative: (1/r)f' ≈ (1/r_i) * (f_{i+1} - f_{i-1})/(2dr)
#             if r_i > 0:
#                 # Coefficient for f_{i-1}
#                 rows.append(i); cols.append(i-1); data.append(-1.0 / (2 * r_i * dr))
#                 # Coefficient for f_{i+1}
#                 rows.append(i); cols.append(i+1); data.append(1.0 / (2 * r_i * dr))
        
#         # Boundary at r=0 (i=0): use symmetry (f' = 0)
#         # For r=0: ∇²f ≈ 2(f_1 - f_0)/dr²
#         rows.append(0); cols.append(0); data.append(-2.0 / dr**2)
#         rows.append(0); cols.append(1); data.append(2.0 / dr**2)
        
#         # Boundary at r=R_max (i=N-1): Neumann boundary (f' = 0)
#         # For r=R_max: ∇²f ≈ 2(f_{N-2} - f_{N-1})/dr²
#         rows.append(N-1); cols.append(N-2); data.append(2.0 / dr**2)
#         rows.append(N-1); cols.append(N-1); data.append(-2.0 / dr**2)
        
#         # Create sparse matrix
#         D2 = sp.coo_matrix((data, (rows, cols)), shape=(N, N), dtype=np.complex128).tocsr()
        
#         return D2
    
#     def background_fields(self, t):
#         cos_wt = np.cos(self.omega * t)
#         sin_wt = np.sin(self.omega * t)
        
#         phi = np.zeros((self.n_fields, self.Nr), dtype=np.complex128)
#         phi[0, :] = self.f0 * np.cos(self.alpha) + self.f1 * cos_wt
#         phi[2, :] = -self.f0 * np.sin(self.alpha) - self.f1 * sin_wt
#         phi[1, :] = -self.f0 * np.cos(self.alpha) + self.f1 * cos_wt
#         phi[3, :] = self.f0 * np.sin(self.alpha) - self.f1 * sin_wt
        
#         return phi
    
#     def compute_potential_jacobian(self, t):
#         phi = self.background_fields(t)
#         psi1, psi2, psi5, psi6 = phi[0, :], phi[1, :], phi[2, :], phi[3, :]
        
#         psi1_sq, psi2_sq, psi5_sq, psi6_sq = psi1**2, psi2**2, psi5**2, psi6**2
#         sum12 = psi1_sq + psi5_sq
#         sum34 = psi2_sq + psi6_sq
#         cross = psi1 * psi2 + psi5 * psi6
        
#         J = np.zeros((self.n_fields, self.n_fields, self.Nr), dtype=np.complex128)
#         mu, lam, lam4 = self.mu, self.lam, self.lam4
        
#         for k in range(self.Nr):
#             p1, p2, p5, p6 = psi1[k], psi2[k], psi5[k], psi6[k]
#             s12, s34, c = sum12[k], sum34[k], cross[k]
            
#             J[0, 0, k] = mu - 2*lam*(3*p1**2 + p5**2 + s34) - 2*lam4*p2**2
#             J[0, 1, k] = -4*lam*p1*p2 - 2*lam4*(c + p1*p2)
#             J[0, 2, k] = -4*lam*p1*p5 - 2*lam4*p2*p6
#             J[0, 3, k] = -4*lam*p1*p6 - 2*lam4*p2*p5
            
#             J[1, 0, k] = -4*lam*p1*p2 - 2*lam4*(c + p1*p2)
#             J[1, 1, k] = mu - 2*lam*(3*p2**2 + p6**2 + s12) - 2*lam4*p1**2
#             J[1, 2, k] = -4*lam*p2*p5 - 2*lam4*p1*p6
#             J[1, 3, k] = -4*lam*p2*p6 - 2*lam4*p1*p5
            
#             J[2, 0, k] = -4*lam*p1*p5 - 2*lam4*p2*p6
#             J[2, 1, k] = -4*lam*p2*p5 - 2*lam4*p1*p6
#             J[2, 2, k] = mu - 2*lam*(3*p5**2 + p1**2 + s34) - 2*lam4*p6**2
#             J[2, 3, k] = -4*lam*p5*p6 - 2*lam4*(c + p5*p6)
            
#             J[3, 0, k] = -4*lam*p1*p6 - 2*lam4*p2*p5
#             J[3, 1, k] = -4*lam*p2*p6 - 2*lam4*p1*p5
#             J[3, 2, k] = -4*lam*p5*p6 - 2*lam4*(c + p5*p6)
#             J[3, 3, k] = mu - 2*lam*(3*p6**2 + p2**2 + s12) - 2*lam4*p5**2
        
#         return J

# # ================================
# # 2. Floquet Analyzer
# # ================================
# class FloquetAnalyzer:
#     def __init__(self, system):
#         self.system = system
#         self.Nr = system.Nr
#         self.n_fields = system.n_fields
#         self.total_dof = system.total_dof
    
#     def build_linear_operator(self, t):
#         # Use 2D radial Laplacian for each field
#         D_list = [self.system.D2_radial for _ in range(self.n_fields)]
#         D_block = sp.block_diag(D_list, format='csr', dtype=np.complex128)
        
#         J_tx = self.system.compute_potential_jacobian(t)
#         rows, cols, data = [], [], []
        
#         for k in range(self.Nr):
#             for i in range(self.n_fields):
#                 for j in range(self.n_fields):
#                     val = J_tx[i, j, k]
#                     if abs(val) > 1e-12:
#                         row_idx = i * self.Nr + k
#                         col_idx = j * self.Nr + k
#                         rows.append(row_idx)
#                         cols.append(col_idx)
#                         data.append(val)
        
#         J_block = sp.coo_matrix((data, (rows, cols)), 
#                                shape=(self.total_dof, self.total_dof),
#                                dtype=np.complex128).tocsr()
        
#         return D_block + J_block
    
#     def compute_monodromy_matrix_smart(self, n_time_samples=10, num_columns=100):
#         T = self.system.T
#         dim = 2 * self.total_dof
        
#         times = np.linspace(0, T, n_time_samples)
#         L_operators = [self.build_linear_operator(t) for t in times]
        
#         M = np.eye(dim, dtype=np.complex128)
#         max_columns = min(num_columns, dim)
        
#         for col in range(max_columns):
#             y0 = M[:, col]
            
#             def ode_rhs(t, y):
#                 idx = min(int(t / T * (n_time_samples - 1)), n_time_samples - 2)
#                 t1, t2 = times[idx], times[idx + 1]
#                 alpha = (t - t1) / (t2 - t1) if t2 > t1 else 0
#                 L = L_operators[idx] if alpha < 0.5 else L_operators[idx + 1]
                
#                 q = y[:self.total_dof]
#                 p = y[self.total_dof:]
#                 return np.concatenate([p, L.dot(q)])
            
#             try:
#                 sol = solve_ivp(ode_rhs, [0, T], y0,
#                               method='RK45',
#                               rtol=1e-8,
#                               atol=1e-10,
#                               max_step=T/500)
                
#                 M[:, col] = sol.y[:, -1]
#             except:
#                 # Fallback integration method
#                 dt = T / 200
#                 y = y0.copy()
#                 for step in range(200):
#                     t_step = step * dt
#                     idx = min(int(t_step / T * (n_time_samples - 1)), n_time_samples - 2)
#                     L = L_operators[idx]
#                     q = y[:self.total_dof]
#                     p = y[self.total_dof:]
#                     dq_dt = p
#                     dp_dt = L.dot(q)
#                     y += dt * np.concatenate([dq_dt, dp_dt])
#                 M[:, col] = y
        
#         return M
    
#     def compute_stability(self, num_modes=100):
#         M = self.compute_monodromy_matrix_smart(num_columns=min(150, 2*self.total_dof//10))
        
#         M_sparse = sp.csr_matrix(M)
#         max_eigenvalues = min(num_modes, M.shape[0]//2)
        
#         try:
#             eigenvalues, _ = spla.eigs(M_sparse, k=max_eigenvalues,
#                                       which='LM',
#                                       return_eigenvectors=True,
#                                       maxiter=10000,
#                                       tol=1e-8)
#         except:
#             print("Warning: Could not compute all requested eigenvalues. Computing fewer...")
#             eigenvalues, _ = spla.eigs(M_sparse, k=min(50, max_eigenvalues),
#                                       which='LM',
#                                       return_eigenvectors=True,
#                                       maxiter=10000,
#                                       tol=1e-6)
        
#         return eigenvalues

# # ================================
# # 3. Plot Functions (Paper Format)
# # ================================
# def plot_background_fields(system):
#     """Plot background fields f0 and f1 in 2D polar coordinates"""
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
#     # Plot f0
#     ax1.plot(system.r, system.f0, 'b-', linewidth=2)
#     ax1.set_xlabel('Radial coordinate r', fontsize=12)
#     ax1.set_ylabel('Field amplitude f₀(r)', fontsize=12)
#     ax1.set_title('Background field f₀(r) in 2D', fontsize=14)
#     ax1.grid(True, alpha=0.3)
#     ax1.set_xlim([0, system.R_max/2])
    
#     # Plot f1
#     ax2.plot(system.r, system.f1, 'r-', linewidth=2)
#     ax2.set_xlabel('Radial coordinate r', fontsize=12)
#     ax2.set_ylabel('Field amplitude f₁(r)', fontsize=12)
#     ax2.set_title('Background field f₁(r) in 2D', fontsize=14)
#     ax2.grid(True, alpha=0.3)
#     ax2.set_xlim([0, system.R_max/2])
    
#     plt.tight_layout()
#     plt.savefig('background_fields_2d.png', dpi=150, bbox_inches='tight')
#     plt.show()

# def create_paper_format_plots(system, eigenvalues):
#     """Create publication-quality plots for 2D polar coordinates - only top row"""
    
#     # Create figure with only 2 subplots (top row)
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
#     # Plot 1: Unit circle plot
#     theta = np.linspace(0, 2*np.pi, 200)
#     ax1.plot(np.cos(theta), np.sin(theta), 'k-', linewidth=1.0, alpha=0.7, label='Unit circle')
    
#     # Separate stable and unstable eigenvalues
#     moduli = np.abs(eigenvalues)
#     stable_mask = moduli <= 0.99999
#     one_mask =  (moduli > 0.99999) & (moduli < 1.00001)
#     unstable_mask = moduli > 1.00001
    
#     stable_eigenvalues = eigenvalues[stable_mask]
#     one_eigenvalues = eigenvalues[one_mask]
#     unstable_eigenvalues = eigenvalues[unstable_mask]
    
#     ax1.scatter(np.real(stable_eigenvalues), np.imag(stable_eigenvalues), 
#                 c='blue', alpha=0.7, edgecolors='none', s=50, label=f'Stable(|ρ|<1): {len(stable_eigenvalues)}')


#     ax1.scatter(np.real(one_eigenvalues), np.imag(one_eigenvalues), 
#                 c='red', alpha=0.9, edgecolors='none', s=50, label=f'stable(|ρ|=1): {len(one_eigenvalues)}')
    

#     ax1.scatter(np.real(unstable_eigenvalues), np.imag(unstable_eigenvalues), 
#                 c='red', alpha=0.9, edgecolors='none', s=70, marker='^', label=f'Unstable(|ρ|>1): {len(unstable_eigenvalues)}')
    
    
#     ax1.set_xlabel('Re(ρ)', fontsize=14)
#     ax1.set_ylabel('Im(ρ)', fontsize=14)
#     ax1.set_title('Floquet Multipliers (2D Polar Coordinates)', fontsize=16)
#     ax1.legend(loc='upper right', fontsize=10)
#     ax1.grid(True, alpha=0.3, linestyle='--')
#     ax1.axis('equal')
#     ax1.set_xlim([-1.5, 1.5])
#     ax1.set_ylim([-1.5, 1.5])
    
#     # Plot 2: Histogram of eigenvalue moduli
#     bins = np.linspace(0, max(2, np.max(moduli)*1.1), 30)
#     n, bins, patches = ax2.hist(moduli, bins=bins, alpha=0.7, color='steelblue', 
#                                edgecolor='black', linewidth=0.5)
    
#     # Color unstable region
#     for i, patch in enumerate(patches):
#         if bins[i] > 1.0:
#             patch.set_facecolor('lightcoral')
#             patch.set_alpha(0.8)
    
#     # Add stability threshold line
#     # ax2.axvline(x=1.0, color='red', linestyle='--', linewidth=2, alpha=0.8, label='Stability threshold')
    
#     ax2.set_xlabel('|ρ|', fontsize=14)
#     ax2.set_ylabel('Frequency', fontsize=14)
#     ax2.set_title('Distribution of Floquet Multiplier Moduli', fontsize=16)
#     # ax2.legend(loc='upper right', fontsize=10)
#     ax2.grid(True, alpha=0.3, linestyle='--')
    
#     # Add statistics text to histogram plot
#     stats_text = f'Total eigenvalues: {len(eigenvalues)}\n'
#     stats_text += f'Stable (|ρ| < 1): {len(stable_eigenvalues)}\n'
#     stats_text += f'Stable (|ρ| = 1): {len(one_eigenvalues)}\n'
#     stats_text += f'Unstable (|ρ| > 1): {len(unstable_eigenvalues)}'
    
#     ax2.text(0.65, 0.95, stats_text, transform=ax2.transAxes,
#             fontsize=10, verticalalignment='top',
#             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
    
#     plt.tight_layout()
#     # plt.savefig('floquet_analysis_2d_simple.pdf', dpi=300)
#     plt.savefig('floquet_analysis_2d_simple.png', dpi=300)
    
#     print(f"\nAnalysis completed:")
#     print(f"  Total eigenvalues analyzed: {len(eigenvalues)}")
#     print(f"  Stable eigenvalues (|ρ| < 1): {len(stable_eigenvalues)}")
#     print(f"  Stable eigenvalues (|ρ| = 1): {len(one_eigenvalues)}")
#     print(f"  Unstable eigenvalues (|ρ| > 1): {len(unstable_eigenvalues)}")
#     if len(unstable_eigenvalues) > 0:
#         print(f"  Maximum modulus: {np.max(moduli):.3f}")
    
#     plt.show()
# # ================================
# # 4. Main Function
# # ================================
# def main():
#     print("=" * 60)
#     print("Floquet Stability Analysis in 2D Polar Coordinates")
#     print(f"Parameters: μ=0.5, λ=1.78, λ₄=-2.56, ω=0.45")
#     print("=" * 60)
    
#     # Initialize system
#     print("\nInitializing system in 2D polar coordinates...")
#     system = OscillonSystem()
    
#     # Plot background fields
#     print("\nPlotting background fields...")
#     plot_background_fields(system)
    
#     # Create analyzer and compute eigenvalues
#     print("\nCreating Floquet analyzer...")
#     analyzer = FloquetAnalyzer(system)
    
#     print("Calculating eigenvalues...")
#     eigenvalues = analyzer.compute_stability(num_modes=100)
    
#     # Create paper-format plots
#     print("\nGenerating plots...")
#     create_paper_format_plots(system, eigenvalues)

# # ================================
# # 5. Run Program
# # ================================
# if __name__ == "__main__":
#     main()











# import numpy as np
# import scipy.sparse as sp
# import scipy.sparse.linalg as spla
# import matplotlib.pyplot as plt
# from scipy.integrate import solve_ivp
# from scipy.optimize import fsolve
# import warnings
# warnings.filterwarnings('ignore')

# # ================================
# # 1. System Parameters and Background Solution
# # ================================
# class OscillonSystem:
#     def __init__(self):
#         # Physical parameters
#         self.mu = 0.5
#         self.lam = 1.78
#         self.lam4 = -2.56
#         self.omega = 0.45
        
#         # Spatial parameters in 2D polar coordinates
#         self.R_max = 32.0  # Maximum radius
#         self.Nr = 51      # Number of radial points
#         self.r = np.linspace(0, self.R_max, self.Nr)
#         self.dr = self.r[1] - self.r[0]
        
#         # Avoid division by zero at r=0
#         self.r_safe = self.r.copy()
#         self.r_safe[0] = 1e-10
        
#         # Time period
#         self.T = 2 * np.pi / self.omega
        
#         # Number of fields
#         self.n_fields = 4
        
#         # Total degrees of freedom
#         self.total_dof = self.n_fields * self.Nr
        
#         # Solve for background fields
#         self.f0, self.f1 = self._solve_background_fields()
#         self.alpha = 0.0
        
#         # Construct 2D radial Laplacian matrix
#         self.D2_radial = self._build_radial_laplacian()
    
#     def _solve_background_fields(self):
#         """Solve for background fields f0 and f1 in 2D polar coordinates"""
#         print("Solving background fields in 2D polar coordinates...")
        
#         N = self.Nr
#         dr = self.dr
#         r = self.r
#         r_safe = self.r_safe
        
#         def equations_2d(y):
#             f0 = y[:N]
#             f1 = y[N:]
#             res0 = np.zeros(N)
#             res1 = np.zeros(N)
            
#             # 2D radial Laplacian: ∇²f = f'' + (1/r)f'
#             for i in range(N):
#                 # Boundary conditions at r=0
#                 if i == 0:
#                     # At r=0: f'(0) = 0, so ∇²f = 2 * (f1 - f0) / dr²
#                     laplacian_f0 = 2 * (f0[1] - f0[0]) / dr**2
#                     laplacian_f1 = 2 * (f1[1] - f1[0]) / dr**2
#                 # Boundary at r=R_max
#                 elif i == N-1:
#                     # Neumann boundary: f'(R_max) = 0
#                     laplacian_f0 = 2 * (f0[N-2] - f0[N-1]) / dr**2
#                     laplacian_f1 = 2 * (f1[N-2] - f1[N-1]) / dr**2
#                 else:
#                     # Interior points
#                     # f'' ≈ (f_{i-1} - 2f_i + f_{i+1}) / dr²
#                     second_deriv_f0 = (f0[i-1] - 2*f0[i] + f0[i+1]) / dr**2
#                     second_deriv_f1 = (f1[i-1] - 2*f1[i] + f1[i+1]) / dr**2
                    
#                     # (1/r)f' ≈ (1/r_i) * (f_{i+1} - f_{i-1})/(2dr)
#                     first_deriv_f0 = (f0[i+1] - f0[i-1]) / (2 * dr)
#                     first_deriv_f1 = (f1[i+1] - f1[i-1]) / (2 * dr)
                    
#                     radial_term_f0 = (1 / r_safe[i]) * first_deriv_f0
#                     radial_term_f1 = (1 / r_safe[i]) * first_deriv_f1
                    
#                     laplacian_f0 = second_deriv_f0 + radial_term_f0
#                     laplacian_f1 = second_deriv_f1 + radial_term_f1
                
#                 # Nonlinear terms
#                 nonlinear_f0 = self.mu * f0[i] - 2*self.lam * f0[i] * (2 * f0[i]**2 + 2 * f1[i]**2) + 2 * self.lam4 * f0[i] * (-f0[i]**2 + f1[i]**2)
#                 nonlinear_f1 = self.omega**2 * f1[i] + self.mu * f1[i] - 2*self.lam * f1[i] * (2 * f0[i]**2 + 2 * f1[i]**2) - 2 * self.lam4 * f1[i] * (-f0[i]**2 + f1[i]**2)
                
#                 res0[i] = laplacian_f0 + nonlinear_f0
#                 res1[i] = laplacian_f1 + nonlinear_f1
            
#             return np.concatenate([res0, res1])
        
#         # Initial guess - Gaussian profiles
        
#         f0_initial = (1+np.tanh((self.r-3)/2))*0.25
#         f1_initial = (1-np.tanh((self.r-3)/2))*0.29
#         y_initial = np.concatenate([f0_initial, f1_initial])
        
#         print("Solving nonlinear equations for f₀ and f₁...")
#         try:
#             solution = fsolve(equations_2d, y_initial, xtol=1e-10, maxfev=50000)
#             print("Background solution found!")
#         except Exception as e:
#             print(f"Warning: fsolve failed with error: {e}")
#             print("Using initial guess as solution")
#             solution = y_initial
        
#         # Extract solutions
#         f0_solution = solution[:self.Nr]
#         f1_solution = solution[self.Nr:]
        

        
#         return f0_solution, f1_solution
    
#     def _build_radial_laplacian(self):
#         """Build 2D radial Laplacian matrix: ∇² = d²/dr² + (1/r) d/dr"""
#         print("Building 2D radial Laplacian matrix...")
        
#         N = self.Nr
#         dr = self.dr
        
#         # Initialize sparse matrix components
#         rows, cols, data = [], [], []
        
#         # Interior points (i=1 to N-2)
#         for i in range(1, N-1):
#             r_i = self.r_safe[i]
            
#             # Second derivative: f'' ≈ (f_{i-1} - 2f_i + f_{i+1})/dr²
#             # Coefficient for f_{i-1}
#             rows.append(i); cols.append(i-1); data.append(1.0 / dr**2)
#             # Coefficient for f_i
#             rows.append(i); cols.append(i); data.append(-2.0 / dr**2)
#             # Coefficient for f_{i+1}
#             rows.append(i); cols.append(i+1); data.append(1.0 / dr**2)
            
#             # First derivative: (1/r)f' ≈ (1/r_i) * (f_{i+1} - f_{i-1})/(2dr)
#             if r_i > 0:
#                 # Coefficient for f_{i-1}
#                 rows.append(i); cols.append(i-1); data.append(-1.0 / (2 * r_i * dr))
#                 # Coefficient for f_{i+1}
#                 rows.append(i); cols.append(i+1); data.append(1.0 / (2 * r_i * dr))
        
#         # Boundary at r=0 (i=0): use symmetry (f' = 0)
#         # For r=0: ∇²f ≈ 2(f_1 - f_0)/dr²
#         rows.append(0); cols.append(0); data.append(-2.0 / dr**2)
#         rows.append(0); cols.append(1); data.append(2.0 / dr**2)
        
#         # Boundary at r=R_max (i=N-1): Neumann boundary (f' = 0)
#         # For r=R_max: ∇²f ≈ 2(f_{N-2} - f_{N-1})/dr²
#         rows.append(N-1); cols.append(N-2); data.append(2.0 / dr**2)
#         rows.append(N-1); cols.append(N-1); data.append(-2.0 / dr**2)
        
#         # Create sparse matrix
#         D2 = sp.coo_matrix((data, (rows, cols)), shape=(N, N), dtype=np.complex128).tocsr()
        
#         return D2
    
#     def background_fields(self, t):
#         cos_wt = np.cos(self.omega * t)
#         sin_wt = np.sin(self.omega * t)
        
#         phi = np.zeros((self.n_fields, self.Nr), dtype=np.complex128)
#         phi[0, :] = self.f0 * np.cos(self.alpha) + self.f1 * cos_wt
#         phi[2, :] = -self.f0 * np.sin(self.alpha) - self.f1 * sin_wt
#         phi[1, :] = -self.f0 * np.cos(self.alpha) + self.f1 * cos_wt
#         phi[3, :] = self.f0 * np.sin(self.alpha) - self.f1 * sin_wt
        
#         return phi
    
#     def compute_potential_jacobian(self, t):
#         phi = self.background_fields(t)
#         psi1, psi2, psi5, psi6 = phi[0, :], phi[1, :], phi[2, :], phi[3, :]
        
#         psi1_sq, psi2_sq, psi5_sq, psi6_sq = psi1**2, psi2**2, psi5**2, psi6**2
#         sum12 = psi1_sq + psi5_sq
#         sum34 = psi2_sq + psi6_sq
#         cross = psi1 * psi2 + psi5 * psi6
        
#         J = np.zeros((self.n_fields, self.n_fields, self.Nr), dtype=np.complex128)
#         mu, lam, lam4 = self.mu, self.lam, self.lam4
        
#         for k in range(self.Nr):
#             p1, p2, p5, p6 = psi1[k], psi2[k], psi5[k], psi6[k]
#             s12, s34, c = sum12[k], sum34[k], cross[k]
            
#             J[0, 0, k] = mu - 2*lam*(3*p1**2 + p5**2 + s34) - 2*lam4*p2**2
#             J[0, 1, k] = -4*lam*p1*p2 - 2*lam4*(c + p1*p2)
#             J[0, 2, k] = -4*lam*p1*p5 - 2*lam4*p2*p6
#             J[0, 3, k] = -4*lam*p1*p6 - 2*lam4*p2*p5
            
#             J[1, 0, k] = -4*lam*p1*p2 - 2*lam4*(c + p1*p2)
#             J[1, 1, k] = mu - 2*lam*(3*p2**2 + p6**2 + s12) - 2*lam4*p1**2
#             J[1, 2, k] = -4*lam*p2*p5 - 2*lam4*p1*p6
#             J[1, 3, k] = -4*lam*p2*p6 - 2*lam4*p1*p5
            
#             J[2, 0, k] = -4*lam*p1*p5 - 2*lam4*p2*p6
#             J[2, 1, k] = -4*lam*p2*p5 - 2*lam4*p1*p6
#             J[2, 2, k] = mu - 2*lam*(3*p5**2 + p1**2 + s34) - 2*lam4*p6**2
#             J[2, 3, k] = -4*lam*p5*p6 - 2*lam4*(c + p5*p6)
            
#             J[3, 0, k] = -4*lam*p1*p6 - 2*lam4*p2*p5
#             J[3, 1, k] = -4*lam*p2*p6 - 2*lam4*p1*p5
#             J[3, 2, k] = -4*lam*p5*p6 - 2*lam4*(c + p5*p6)
#             J[3, 3, k] = mu - 2*lam*(3*p6**2 + p2**2 + s12) - 2*lam4*p5**2
        
#         return J

# # ================================
# # 2. Floquet Analyzer
# # ================================
# class FloquetAnalyzer:
#     def __init__(self, system):
#         self.system = system
#         self.Nr = system.Nr
#         self.n_fields = system.n_fields
#         self.total_dof = system.total_dof
    
#     def build_linear_operator(self, t):
#         """Build full linear operator with Sponge Layer ABC: d/dt(y) = M * y
#            where y = [delta, dot_delta].
#            M = [[0, I], [D2+J, -Sponge]]
#            The 'Sponge' term introduces the delta_dot dependence (absorption)."""
        
#         # --- 1. Construct the original Force Matrix (Laplacian + Jacobian) ---
#         # Use 2D radial Laplacian for each field
#         D_list = [self.system.D2_radial for _ in range(self.n_fields)]
#         D_block = sp.block_diag(D_list, format='csr', dtype=np.complex128)
        
#         J_tx = self.system.compute_potential_jacobian(t)
#         rows, cols, data = [], [], []
        
#         for k in range(self.Nr):
#             for i in range(self.n_fields):
#                 for j in range(self.n_fields):
#                     val = J_tx[i, j, k]
#                     if abs(val) > 1e-13:
#                         row_idx = i * self.Nr + k
#                         col_idx = j * self.Nr + k
#                         rows.append(row_idx)
#                         cols.append(col_idx)
#                         data.append(val)
        
#         J_block = sp.coo_matrix((data, (rows, cols)), 
#                                shape=(self.total_dof, self.total_dof),
#                                dtype=np.complex128).tocsr()
        
#         L_force = D_block + J_block
        
#         # --- 2. Construct the Sponge Layer (Absorbing Boundary Condition) ---
#         # This adds the term -Gamma(r) * dot_delta
#         r = self.system.r
#         R_max = self.system.R_max
        
#         # Sponge parameters
#         absorb_width = 3.0
#         absorb_start = R_max - absorb_width
#         sponge_strength = 20.0  # Adjust strength as needed
        
#         # Define sponge profile (0 in bulk, quadratic/tanh rise near boundary)
#         sponge_profile = np.zeros_like(r)
#         mask = r > absorb_start
#         if np.any(mask):
#             # Quadratic ramp
#             sponge_profile[mask] = sponge_strength * ((r[mask] - absorb_start) / absorb_width)**2
        
#         # Create diagonal sponge matrix for one field
#         S_diag = sp.diags(sponge_profile, format='csr', dtype=np.complex128)
        
#         # Block diagonal sponge for all fields
#         S_block = sp.block_diag([S_diag] * self.n_fields, format='csr')
        
#         # --- 3. Assemble the Full 2N x 2N Phase Space Matrix ---
#         # Matrix form:
#         # [ 0       I     ]
#         # [ L_force -S_block ]
        
#         dim = self.total_dof
#         Zero = sp.csr_matrix((dim, dim), dtype=np.complex128)
#         Eye = sp.eye(dim, format='csr', dtype=np.complex128)
        
#         # Use bmat to assemble the sparse blocks efficiently
#         Full_L = sp.bmat([
#             [Zero,    Eye],
#             [L_force, -S_block]
#         ], format='csr')
        
#         return Full_L
    
#     def compute_full_monodromy_matrix_high_precision(self, n_time_samples=50):
#         """Compute full monodromy matrix with high precision using the full operator"""
#         T = self.system.T
#         dim = 2 * self.total_dof
        
#         print(f"Computing full monodromy matrix with high precision (with Sponge Layer)...")
#         print(f"Dimension: {dim} × {dim}")
#         print(f"Using {n_time_samples} time samples")
        
#         # Precompute linear operators at different times with higher sampling
#         times = np.linspace(0, T, n_time_samples, endpoint=False)
#         L_operators = []
        
#         print("Precomputing linear operators at time samples...")
#         for i, t in enumerate(times):
#             if i % 1000 == 0:
#                 print(f"  Time sample {i+1}/{n_time_samples}: t = {t:.4f}")
#             L_operators.append(self.build_linear_operator(t))
        
#         # Initialize monodromy matrix as identity
#         M = np.eye(dim, dtype=np.complex128)
        
#         # Compute all columns with high precision
#         print("Computing monodromy matrix columns...")
#         for col in range(dim):
#             if col % 100 == 0:
#                 print(f"  Column {col+1}/{dim}")
            
#             y0 = M[:, col]
            
#             def ode_rhs(t, y):
#                 # Find nearest time sample with linear interpolation
#                 idx = int(t / T * n_time_samples) % n_time_samples
#                 idx_next = (idx + 1) % n_time_samples
                
#                 t1 = times[idx]
#                 t2 = times[idx_next] if idx_next != 0 else T
                
#                 # Linear interpolation between operators
#                 if t2 > t1:
#                     alpha = (t - t1) / (t2 - t1)
#                 else:
#                     alpha = (t - t1) / (T - t1 + t2)
                
#                 # Interpolate between the two nearest operators
#                 # Note: L is now the full 2N x 2N matrix including sponge terms
#                 L1 = L_operators[idx]
#                 L2 = L_operators[idx_next] if idx_next < n_time_samples else L_operators[0]
#                 L = (1.0 - alpha) * L1 + alpha * L2
                
#                 # The equation is simply dy/dt = L * y
#                 return L.dot(y)
            
#             try:
#                 # Use high-order integration method with very tight tolerances
#                 sol = solve_ivp(ode_rhs, [0, T], y0,
#                               method='DOP853',
#                               rtol=1e-9,
#                               atol=1e-11,
#                               max_step=T/200,
#                               dense_output=False)
                
#                 M[:, col] = sol.y[:, -1]
                
#             except Exception as e:
#                 print(f"  Warning: Integration failed for column {col}: {e}")
#                 print(f"  Using identity column as fallback")
        
#         print("Monodromy matrix computation completed!")
#         return M
    
#     def compute_all_eigenvalues_high_precision(self):
#         """Compute all eigenvalues with high precision"""
#         print("\nComputing full monodromy matrix with high precision...")
#         M = self.compute_full_monodromy_matrix_high_precision(n_time_samples=50000)
        
#         print("\nComputing all eigenvalues and eigenvectors with high precision...")
#         print("This may take a long time, please be patient...")
        
#         # Convert to dense matrix for eigenvalue computation
#         M_dense = M.toarray() if sp.issparse(M) else M
        
#         # Use numpy's high-precision eigenvalue solver
#         eigenvalues, eigenvectors = np.linalg.eig(M_dense)
        
#         # Sort eigenvalues and eigenvectors by magnitude for better analysis
#         sorted_indices = np.argsort(-np.abs(eigenvalues))
#         eigenvalues = eigenvalues[sorted_indices]
#         eigenvectors = eigenvectors[:, sorted_indices]
        
#         print(f"Computation completed! Found {len(eigenvalues)} eigenvalues and eigenvectors")

#         return eigenvalues, eigenvectors, M 
    
#     def analyze_stability_high_precision(self):
#         """Complete stability analysis with high precision"""
#         eigenvalues, eigenvectors, M = self.compute_all_eigenvalues_high_precision()
        
#         # Compute Floquet exponents
#         T = self.system.T
#         exponents = np.log(eigenvalues) / T
        
#         # Analyze results with high precision
#         moduli = np.abs(eigenvalues)
        
#         # Use tighter thresholds for stability classification
#         stable_mask = moduli <= 0.99999
#         marginal_mask = (moduli > 0.99999) & (moduli < 1.00001)
#         unstable_mask = moduli >= 1.00001
        
#         results = {
#             'multipliers': eigenvalues,
#             'exponents': exponents,
#             'eigenvectors': eigenvectors,
#             'moduli': moduli,
#             'num_total': len(eigenvalues),
#             'num_stable': np.sum(stable_mask),
#             'num_marginal': np.sum(marginal_mask),
#             'num_unstable': np.sum(unstable_mask),
#             'stable_indices': np.where(stable_mask)[0],
#             'marginal_indices': np.where(marginal_mask)[0],
#             'unstable_indices': np.where(unstable_mask)[0],
#             'monodromy_matrix': M
#         }
        
#         if np.sum(unstable_mask) > 0:
#             unstable_growth_rates = np.real(exponents[unstable_mask])
#             results['max_growth_rate'] = np.max(unstable_growth_rates)
#             results['avg_growth_rate'] = np.mean(unstable_growth_rates)
#             results['max_modulus'] = np.max(moduli[unstable_mask])
            
#             # Sort unstable modes by growth rate
#             unstable_indices = results['unstable_indices']
#             sorted_unstable_idx = unstable_indices[np.argsort(-unstable_growth_rates)]
#             results['sorted_unstable_indices'] = sorted_unstable_idx
            
#             # Also store the most unstable mode index
#             results['most_unstable_index'] = sorted_unstable_idx[0]
#         else:
#             results['max_growth_rate'] = 0.0
#             results['avg_growth_rate'] = 0.0
#             results['max_modulus'] = np.max(moduli)
#             results['sorted_unstable_indices'] = np.array([])
#             results['most_unstable_index'] = None
        
#         return results


# # ================================
# # 3. Plot Functions (Paper Format)
# # ================================
# def plot_background_fields(system):
#     """Plot background fields f0 and f1 in 2D polar coordinates"""
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
#     # Plot f0
#     ax1.plot(system.r, system.f0, 'b-', linewidth=2)
#     ax1.set_xlabel('Radial coordinate r', fontsize=12)
#     ax1.set_ylabel('Field amplitude f₀(r)', fontsize=12)
#     ax1.set_title('Background field f₀(r) in 2D', fontsize=14)
#     ax1.grid(True, alpha=0.3)
#     ax1.set_xlim([0, system.R_max/2])
    
#     # Plot f1
#     ax2.plot(system.r, system.f1, 'r-', linewidth=2)
#     ax2.set_xlabel('Radial coordinate r', fontsize=12)
#     ax2.set_ylabel('Field amplitude f₁(r)', fontsize=12)
#     ax2.set_title('Background field f₁(r) in 2D', fontsize=14)
#     ax2.grid(True, alpha=0.3)
#     ax2.set_xlim([0, system.R_max/2])
    
#     plt.tight_layout()
#     plt.savefig('background_fields_2d.png', dpi=150, bbox_inches='tight')
#     plt.show()

# def create_paper_format_plots(system, results):
#     """Create publication-quality plots for 2D polar coordinates - only top row"""
    

#     eigenvalues=results['multipliers']
#     # Create figure with only 2 subplots (top row)
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
#     # Plot 1: Unit circle plot
#     theta = np.linspace(0, 2*np.pi, 200)
#     ax1.plot(np.cos(theta), np.sin(theta), 'k-', linewidth=1.0, alpha=0.7, label='Unit circle')
    
#     # Separate stable and unstable eigenvalues
#     moduli = np.abs(eigenvalues)
#     stable_mask = moduli <= 0.99999
#     one_mask =  (moduli > 0.99999) & (moduli < 1.00001)
#     unstable_mask = moduli > 1.00001
    
#     stable_eigenvalues = eigenvalues[stable_mask]
#     one_eigenvalues = eigenvalues[one_mask]
#     unstable_eigenvalues = eigenvalues[unstable_mask]
    
#     ax1.scatter(np.real(stable_eigenvalues), np.imag(stable_eigenvalues), 
#                 c='blue', alpha=0.7, edgecolors='none', s=50, label=f'Stable(|ρ|<1): {len(stable_eigenvalues)}')


#     ax1.scatter(np.real(one_eigenvalues), np.imag(one_eigenvalues), 
#                 c='red', alpha=0.9, edgecolors='none', s=50, label=f'stable(|ρ|=1): {len(one_eigenvalues)}')
    

#     ax1.scatter(np.real(unstable_eigenvalues), np.imag(unstable_eigenvalues), 
#                 c='red', alpha=0.9, edgecolors='none', s=70, marker='^', label=f'Unstable(|ρ|>1): {len(unstable_eigenvalues)}')
    
    
#     ax1.set_xlabel('Re(ρ)', fontsize=14)
#     ax1.set_ylabel('Im(ρ)', fontsize=14)
#     ax1.set_title('Floquet Multipliers (2D Polar Coordinates)', fontsize=16)
#     ax1.legend(loc='upper right', fontsize=10)
#     ax1.grid(True, alpha=0.3, linestyle='--')
#     ax1.axis('equal')
#     ax1.set_xlim([-1.5, 1.5])
#     ax1.set_ylim([-1.5, 1.5])
    
#     # Plot 2: Histogram of eigenvalue moduli
#     bins = np.linspace(0, max(2, np.max(moduli)*1.1), 30)
#     n, bins, patches = ax2.hist(moduli, bins=bins, alpha=0.7, color='steelblue', 
#                                edgecolor='black', linewidth=0.5)
    
#     # Color unstable region
#     for i, patch in enumerate(patches):
#         if bins[i] > 1.0:
#             patch.set_facecolor('lightcoral')
#             patch.set_alpha(0.8)
    
#     # Add stability threshold line
#     # ax2.axvline(x=1.0, color='red', linestyle='--', linewidth=2, alpha=0.8, label='Stability threshold')
    
#     ax2.set_xlabel('|ρ|', fontsize=14)
#     ax2.set_ylabel('Frequency', fontsize=14)
#     ax2.set_title('Distribution of Floquet Multiplier Moduli', fontsize=16)
#     # ax2.legend(loc='upper right', fontsize=10)
#     ax2.grid(True, alpha=0.3, linestyle='--')
    
#     # Add statistics text to histogram plot
#     stats_text = f'Total eigenvalues: {len(eigenvalues)}\n'
#     stats_text += f'Stable (|ρ| < 1): {len(stable_eigenvalues)}\n'
#     stats_text += f'Stable (|ρ| = 1): {len(one_eigenvalues)}\n'
#     stats_text += f'Unstable (|ρ| > 1): {len(unstable_eigenvalues)}'
    
#     ax2.text(0.65, 0.95, stats_text, transform=ax2.transAxes,
#             fontsize=10, verticalalignment='top',
#             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
    
#     plt.tight_layout()
#     # plt.savefig('floquet_analysis_2d_simple.pdf', dpi=300)
#     plt.savefig('floquet_analysis_2d_simple.png', dpi=300)
    
#     print(f"\nAnalysis completed:")
#     print(f"  Total eigenvalues analyzed: {len(eigenvalues)}")
#     print(f"  Stable eigenvalues (|ρ| < 1): {len(stable_eigenvalues)}")
#     print(f"  Stable eigenvalues (|ρ| = 1): {len(one_eigenvalues)}")
#     print(f"  Unstable eigenvalues (|ρ| > 1): {len(unstable_eigenvalues)}")
#     if len(unstable_eigenvalues) > 0:
#         print(f"  Maximum modulus: {np.max(moduli):.3f}")
    
#     plt.show()
# # ================================
# # 4. Main Function
# # ================================
# def main():
#     print("=" * 60)
#     print("Floquet Stability Analysis in 2D Polar Coordinates")
#     print(f"Parameters: μ=0.5, λ=1.78, λ₄=-2.56, ω=0.45")
#     print("=" * 60)
    
#     # Initialize system
#     print("\nInitializing system in 2D polar coordinates...")
#     system = OscillonSystem()
    
#     # Plot background fields
#     print("\nPlotting background fields...")
#     plot_background_fields(system)
    
#     # Create analyzer and compute eigenvalues
#     print("\nCreating Floquet analyzer...")
#     analyzer = FloquetAnalyzer(system)
    
#     print("Calculating eigenvalues...")
#     # eigenvalues = analyzer.compute_stability(num_modes=100)
#     results = analyzer.analyze_stability_high_precision()
    
#     # Create paper-format plots
#     print("\nGenerating plots...")
#     create_paper_format_plots(system, results)

# # ================================
# # 5. Run Program
# # ================================
# if __name__ == "__main__":
#     main()
















import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import warnings
warnings.filterwarnings('ignore')

# ================================
# 1. System Parameters and Background Solution
# ================================
class OscillonSystem:
    def __init__(self):
        # Physical parameters
        self.mu = 0.5
        self.lam = 1.78
        self.lam4 = -2.56
        self.omega = 0.45
        
        # Spatial parameters in 2D polar coordinates
        self.R_max = 32.0  # Maximum radius
        self.Nr = 129      # Number of radial points
        self.r = np.linspace(0, self.R_max, self.Nr)
        self.dr = self.r[1] - self.r[0]
        
        # Avoid division by zero at r=0
        self.r_safe = self.r.copy()
        self.r_safe[0] = 1e-10
        
        # Time period
        self.T = 2 * np.pi / self.omega
        
        # Number of fields
        self.n_fields = 4
        
        # Total degrees of freedom
        self.total_dof = self.n_fields * self.Nr
        
        # Solve for background fields
        self.f0, self.f1 = self._solve_background_fields()
        self.alpha = 0.0
        
        # Construct 2D radial Laplacian matrix
        self.D2_radial = self._build_radial_laplacian()
    
    def _solve_background_fields(self):
        """Solve for background fields f0 and f1 in 2D polar coordinates"""
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
            
            # 2D radial Laplacian: ∇²f = f'' + (1/r)f'
            for i in range(N):
                # Boundary conditions at r=0
                if i == 0:
                    # At r=0: f'(0) = 0, so ∇²f = 2 * (f1 - f0) / dr²
                    laplacian_f0 = 2 * (f0[1] - f0[0]) / dr**2
                    laplacian_f1 = 2 * (f1[1] - f1[0]) / dr**2
                # Boundary at r=R_max
                elif i == N-1:
                    # Neumann boundary: f'(R_max) = 0
                    laplacian_f0 = 2 * (f0[N-2] - f0[N-1]) / dr**2
                    laplacian_f1 = 2 * (f1[N-2] - f1[N-1]) / dr**2
                else:
                    # Interior points
                    # f'' ≈ (f_{i-1} - 2f_i + f_{i+1}) / dr²
                    second_deriv_f0 = (f0[i-1] - 2*f0[i] + f0[i+1]) / dr**2
                    second_deriv_f1 = (f1[i-1] - 2*f1[i] + f1[i+1]) / dr**2
                    
                    # (1/r)f' ≈ (1/r_i) * (f_{i+1} - f_{i-1})/(2dr)
                    first_deriv_f0 = (f0[i+1] - f0[i-1]) / (2 * dr)
                    first_deriv_f1 = (f1[i+1] - f1[i-1]) / (2 * dr)
                    
                    radial_term_f0 = (1 / r_safe[i]) * first_deriv_f0
                    radial_term_f1 = (1 / r_safe[i]) * first_deriv_f1
                    
                    laplacian_f0 = second_deriv_f0 + radial_term_f0
                    laplacian_f1 = second_deriv_f1 + radial_term_f1
                
                # Nonlinear terms
                nonlinear_f0 = self.mu * f0[i] - 2*self.lam * f0[i] * (2 * f0[i]**2 + 2 * f1[i]**2) + 2 * self.lam4 * f0[i] * (-f0[i]**2 + f1[i]**2)
                nonlinear_f1 = self.omega**2 * f1[i] + self.mu * f1[i] - 2*self.lam * f1[i] * (2 * f0[i]**2 + 2 * f1[i]**2) - 2 * self.lam4 * f1[i] * (-f0[i]**2 + f1[i]**2)
                
                res0[i] = laplacian_f0 + nonlinear_f0
                res1[i] = laplacian_f1 + nonlinear_f1
            
            return np.concatenate([res0, res1])
        
        # Initial guess - Gaussian profiles
        
        f0_initial = (1+np.tanh((self.r-3)/2))*0.25
        f1_initial = (1-np.tanh((self.r-3)/2))*0.29
        y_initial = np.concatenate([f0_initial, f1_initial])
        
        print("Solving nonlinear equations for f₀ and f₁...")
        try:
            solution = fsolve(equations_2d, y_initial, xtol=1e-11, maxfev=5000000)
            print("Background solution found!")
        except Exception as e:
            print(f"Warning: fsolve failed with error: {e}")
            print("Using initial guess as solution")
            solution = y_initial
        
        # Extract solutions
        f0_solution = solution[:self.Nr]
        f1_solution = solution[self.Nr:]
        

        
        return f0_solution, f1_solution
    
    def _build_radial_laplacian(self):
        """Build 2D radial Laplacian matrix: ∇² = d²/dr² + (1/r) d/dr"""
        print("Building 2D radial Laplacian matrix...")
        
        N = self.Nr
        dr = self.dr
        
        # Initialize sparse matrix components
        rows, cols, data = [], [], []
        
        # Interior points (i=1 to N-2)
        for i in range(1, N-1):
            r_i = self.r_safe[i]
            
            # Second derivative: f'' ≈ (f_{i-1} - 2f_i + f_{i+1})/dr²
            # Coefficient for f_{i-1}
            rows.append(i); cols.append(i-1); data.append(1.0 / dr**2)
            # Coefficient for f_i
            rows.append(i); cols.append(i); data.append(-2.0 / dr**2)
            # Coefficient for f_{i+1}
            rows.append(i); cols.append(i+1); data.append(1.0 / dr**2)
            
            # First derivative: (1/r)f' ≈ (1/r_i) * (f_{i+1} - f_{i-1})/(2dr)
            if r_i > 0:
                # Coefficient for f_{i-1}
                rows.append(i); cols.append(i-1); data.append(-1.0 / (2 * r_i * dr))
                # Coefficient for f_{i+1}
                rows.append(i); cols.append(i+1); data.append(1.0 / (2 * r_i * dr))
        
        # Boundary at r=0 (i=0): use symmetry (f' = 0)
        # For r=0: ∇²f ≈ 2(f_1 - f_0)/dr²
        rows.append(0); cols.append(0); data.append(-2.0 / dr**2)
        rows.append(0); cols.append(1); data.append(2.0 / dr**2)
        
        # Boundary at r=R_max (i=N-1): Neumann boundary (f' = 0)
        # For r=R_max: ∇²f ≈ 2(f_{N-2} - f_{N-1})/dr²
        rows.append(N-1); cols.append(N-2); data.append(2.0 / dr**2)
        rows.append(N-1); cols.append(N-1); data.append(-2.0 / dr**2)
        
        # Create sparse matrix
        D2 = sp.coo_matrix((data, (rows, cols)), shape=(N, N), dtype=np.complex128).tocsr()
        
        return D2
    
    def background_fields(self, t):
        cos_wt = np.cos(self.omega * t)
        sin_wt = np.sin(self.omega * t)
        
        phi = np.zeros((self.n_fields, self.Nr), dtype=np.complex128)
        phi[0, :] = self.f0 * np.cos(self.alpha) + self.f1 * cos_wt
        phi[2, :] = -self.f0 * np.sin(self.alpha) - self.f1 * sin_wt
        phi[1, :] = -self.f0 * np.cos(self.alpha) + self.f1 * cos_wt
        phi[3, :] = self.f0 * np.sin(self.alpha) - self.f1 * sin_wt
        
        return phi
    
    def compute_potential_jacobian(self, t):
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
            
            J[0, 0, k] = mu - 2*lam*(3*p1**2 + p5**2 + s34) - 2*lam4*p2**2
            J[0, 1, k] = -4*lam*p1*p2 - 2*lam4*(c + p1*p2)
            J[0, 2, k] = -4*lam*p1*p5 - 2*lam4*p2*p6
            J[0, 3, k] = -4*lam*p1*p6 - 2*lam4*p2*p5
            
            J[1, 0, k] = -4*lam*p1*p2 - 2*lam4*(c + p1*p2)
            J[1, 1, k] = mu - 2*lam*(3*p2**2 + p6**2 + s12) - 2*lam4*p1**2
            J[1, 2, k] = -4*lam*p2*p5 - 2*lam4*p1*p6
            J[1, 3, k] = -4*lam*p2*p6 - 2*lam4*p1*p5
            
            J[2, 0, k] = -4*lam*p1*p5 - 2*lam4*p2*p6
            J[2, 1, k] = -4*lam*p2*p5 - 2*lam4*p1*p6
            J[2, 2, k] = mu - 2*lam*(3*p5**2 + p1**2 + s34) - 2*lam4*p6**2
            J[2, 3, k] = -4*lam*p5*p6 - 2*lam4*(c + p5*p6)
            
            J[3, 0, k] = -4*lam*p1*p6 - 2*lam4*p2*p5
            J[3, 1, k] = -4*lam*p2*p6 - 2*lam4*p1*p5
            J[3, 2, k] = -4*lam*p5*p6 - 2*lam4*(c + p5*p6)
            J[3, 3, k] = mu - 2*lam*(3*p6**2 + p2**2 + s12) - 2*lam4*p5**2
        
        return J

# ================================
# 2. Floquet Analyzer
# ================================
class FloquetAnalyzer:
    def __init__(self, system):
        self.system = system
        self.Nr = system.Nr
        self.n_fields = system.n_fields
        self.total_dof = system.total_dof
    
    def build_linear_operator(self, t):
        """Build full linear operator with Sponge Layer ABC: d/dt(y) = M * y
           where y = [delta, dot_delta].
           M = [[0, I], [D2+J, -Sponge]]
           The 'Sponge' term introduces the delta_dot dependence (absorption)."""
        
        # --- 1. Construct the original Force Matrix (Laplacian + Jacobian) ---
        # Use 2D radial Laplacian for each field
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
        
        # --- 2. Construct the Sponge Layer (Absorbing Boundary Condition) ---
        # This adds the term -Gamma(r) * dot_delta
        r = self.system.r
        R_max = self.system.R_max
        
        # Sponge parameters
        absorb_width = R_max*0.1
        absorb_start = R_max*0.9
        sponge_strength = 5.0  # Adjust strength as needed
        
        # Define sponge profile (0 in bulk, quadratic/tanh rise near boundary)
        sponge_profile = np.zeros_like(r)
        mask = r > absorb_start
        if np.any(mask):
            # Quadratic ramp
            sponge_profile[mask] = sponge_strength * ((r[mask] - absorb_start) / absorb_width)**2
        
        # Create diagonal sponge matrix for one field
        S_diag = sp.diags(sponge_profile, format='csr', dtype=np.complex128)
        
        # Block diagonal sponge for all fields
        S_block = sp.block_diag([S_diag] * self.n_fields, format='csr')
        
        # --- 3. Assemble the Full 2N x 2N Phase Space Matrix ---
        # Matrix form:
        # [ 0       I     ]
        # [ L_force -S_block ]
        
        dim = self.total_dof
        Zero = sp.csr_matrix((dim, dim), dtype=np.complex128)
        Eye = sp.eye(dim, format='csr', dtype=np.complex128)
        
        # Use bmat to assemble the sparse blocks efficiently
        Full_L = sp.bmat([
            [Zero,    Eye],
            [L_force, -S_block]
        ], format='csr')
        
        return Full_L
    
    def compute_full_monodromy_matrix_high_precision(self, n_time_samples=50):
        """Compute full monodromy matrix with high precision using the full operator"""
        T = self.system.T
        dim = 2 * self.total_dof
        
        print(f"Computing full monodromy matrix with high precision (with Sponge Layer)...")
        print(f"Dimension: {dim} × {dim}")
        print(f"Using {n_time_samples} time samples")
        
        # # Precompute linear operators at different times with higher sampling
        # times = np.linspace(0, T, n_time_samples, endpoint=False)
        # L_operators = []
        
        # print("Precomputing linear operators at time samples...")
        # for i, t in enumerate(times):
        #     if i % 1000 == 0:
        #         print(f"  Time sample {i+1}/{n_time_samples}: t = {t:.4f}")
        #     L_operators.append(self.build_linear_operator(t))
        
        # Initialize monodromy matrix as identity
        M = np.eye(dim, dtype=np.complex128)
        
        # Compute all columns with high precision
        print("Computing monodromy matrix columns...")
        for col in range(dim):
            if col % 100 == 0:
                print(f"  Column {col+1}/{dim}")
            
            y0 = M[:, col]
            
            # def ode_rhs(t, y):
            #     # Find nearest time sample with linear interpolation
            #     idx = int(t / T * n_time_samples) % n_time_samples
            #     idx_next = (idx + 1) % n_time_samples
                
            #     t1 = times[idx]
            #     t2 = times[idx_next] if idx_next != 0 else T
                
            #     # Linear interpolation between operators
            #     if t2 > t1:
            #         alpha = (t - t1) / (t2 - t1)
            #     else:
            #         alpha = (t - t1) / (T - t1 + t2)
                
            #     # Interpolate between the two nearest operators
            #     # Note: L is now the full 2N x 2N matrix including sponge terms
            #     L1 = L_operators[idx]
            #     L2 = L_operators[idx_next] if idx_next < n_time_samples else L_operators[0]
            #     L = (1.0 - alpha) * L1 + alpha * L2
                
            #     # The equation is simply dy/dt = L * y
            #     # L=L_operators[idx]
            #     return L.dot(y)
            def ode_rhs(t, y):
                # 1. 直接通过解析公式计算当前时刻的 Jacobian (无插值误差)
                # 2. build_linear_operator 内部调用 background_fields(t)
                L = self.build_linear_operator(t) 
                
                # 3. y 是当前的扰动向量，直接矩阵相乘
                return L.dot(y)
            try:
                # Use high-order integration method with very tight tolerances
                sol = solve_ivp(ode_rhs, [0, T], y0,
                              method='DOP853',
                              rtol=1e-11,
                              atol=1e-13,
                              max_step=T/140,
                              dense_output=False)
                
                M[:, col] = sol.y[:, -1]
                
            except Exception as e:
                print(f"  Warning: Integration failed for column {col}: {e}")
                print(f"  Using identity column as fallback")
        
        print("Monodromy matrix computation completed!")
        return M
    
    def compute_all_eigenvalues_high_precision(self):
        """Compute all eigenvalues with high precision"""
        print("\nComputing full monodromy matrix with high precision...")
        M = self.compute_full_monodromy_matrix_high_precision(n_time_samples=140000)
        
        print("\nComputing all eigenvalues and eigenvectors with high precision...")
        print("This may take a long time, please be patient...")
        
        # Convert to dense matrix for eigenvalue computation
        M_dense = M.toarray() if sp.issparse(M) else M
        
        # Use numpy's high-precision eigenvalue solver
        eigenvalues, eigenvectors = np.linalg.eig(M_dense)
        
        # Sort eigenvalues and eigenvectors by magnitude for better analysis
        sorted_indices = np.argsort(-np.abs(eigenvalues))
        eigenvalues = eigenvalues[sorted_indices]
        eigenvectors = eigenvectors[:, sorted_indices]
        
        print(f"Computation completed! Found {len(eigenvalues)} eigenvalues and eigenvectors")

        return eigenvalues, eigenvectors, M 
    
    def analyze_stability_high_precision(self):
        """Complete stability analysis with high precision"""
        eigenvalues, eigenvectors, M = self.compute_all_eigenvalues_high_precision()
        
        # Compute Floquet exponents
        T = self.system.T
        exponents = np.log(eigenvalues) / T
        
        # Analyze results with high precision
        moduli = np.abs(eigenvalues)
        
        # Use tighter thresholds for stability classification
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
            
            # Also store the most unstable mode index
            results['most_unstable_index'] = sorted_unstable_idx[0]
        else:
            results['max_growth_rate'] = 0.0
            results['avg_growth_rate'] = 0.0
            results['max_modulus'] = np.max(moduli)
            results['sorted_unstable_indices'] = np.array([])
            results['most_unstable_index'] = None
        
        return results


# ================================
# 3. Plot Functions (Paper Format)
# ================================
def plot_background_fields(system):
    """Plot background fields f0 and f1 in 2D polar coordinates"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot f0
    ax1.plot(system.r, system.f0, 'b-', linewidth=2)
    ax1.set_xlabel('Radial coordinate r', fontsize=12)
    ax1.set_ylabel('Field amplitude f₀(r)', fontsize=12)
    ax1.set_title('Background field f₀(r) in 2D', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([0, system.R_max/2])
    
    # Plot f1
    ax2.plot(system.r, system.f1, 'r-', linewidth=2)
    ax2.set_xlabel('Radial coordinate r', fontsize=12)
    ax2.set_ylabel('Field amplitude f₁(r)', fontsize=12)
    ax2.set_title('Background field f₁(r) in 2D', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([0, system.R_max/2])
    
    plt.tight_layout()
    plt.savefig('background_fields_2d_4field.png', dpi=150, bbox_inches='tight')
    plt.show()

def create_paper_format_plots(system, results):
    """Create publication-quality plots for 2D polar coordinates - only top row"""
    

    eigenvalues=results['multipliers']
    # Create figure with only 2 subplots (top row)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: Unit circle plot
    theta = np.linspace(0, 2*np.pi, 200)
    ax1.plot(np.cos(theta), np.sin(theta), 'k-', linewidth=1.0, alpha=0.7, label='Unit circle')
    
    # Separate stable and unstable eigenvalues
    moduli = np.abs(eigenvalues)
    stable_mask = moduli <= 0.99999
    one_mask =  (moduli > 0.99999) & (moduli < 1.00001)
    unstable_mask = moduli > 1.00001
    
    stable_eigenvalues = eigenvalues[stable_mask]
    one_eigenvalues = eigenvalues[one_mask]
    unstable_eigenvalues = eigenvalues[unstable_mask]
    
    ax1.scatter(np.real(stable_eigenvalues), np.imag(stable_eigenvalues), 
                c='blue', alpha=0.7, edgecolors='none', s=50, label=f'Stable(|ρ|<1): {len(stable_eigenvalues)}')


    ax1.scatter(np.real(one_eigenvalues), np.imag(one_eigenvalues), 
                c='red', alpha=0.9, edgecolors='none', s=50, label=f'stable(|ρ|=1): {len(one_eigenvalues)}')
    

    ax1.scatter(np.real(unstable_eigenvalues), np.imag(unstable_eigenvalues), 
                c='red', alpha=0.9, edgecolors='none', s=70, marker='^', label=f'Unstable(|ρ|>1): {len(unstable_eigenvalues)}')
    
    
    ax1.set_xlabel('Re(ρ)', fontsize=14)
    ax1.set_ylabel('Im(ρ)', fontsize=14)
    ax1.set_title('Floquet Multipliers (2D Polar Coordinates)', fontsize=16)
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.axis('equal')
    ax1.set_xlim([-1.5, 1.5])
    ax1.set_ylim([-1.5, 1.5])
    
    # Plot 2: Histogram of eigenvalue moduli
    bins = np.linspace(0, max(2, np.max(moduli)*1.1), 30)
    n, bins, patches = ax2.hist(moduli, bins=bins, alpha=0.7, color='steelblue', 
                               edgecolor='black', linewidth=0.5)
    
    # Color unstable region
    for i, patch in enumerate(patches):
        if bins[i] > 1.0:
            patch.set_facecolor('lightcoral')
            patch.set_alpha(0.8)
    
    # Add stability threshold line
    # ax2.axvline(x=1.0, color='red', linestyle='--', linewidth=2, alpha=0.8, label='Stability threshold')
    
    ax2.set_xlabel('|ρ|', fontsize=14)
    ax2.set_ylabel('Frequency', fontsize=14)
    ax2.set_title('Distribution of Floquet Multiplier Moduli', fontsize=16)
    # ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(True, alpha=0.3, linestyle='--')
    
    # Add statistics text to histogram plot
    stats_text = f'Total eigenvalues: {len(eigenvalues)}\n'
    stats_text += f'Stable (|ρ| < 1): {len(stable_eigenvalues)}\n'
    stats_text += f'Stable (|ρ| = 1): {len(one_eigenvalues)}\n'
    stats_text += f'Unstable (|ρ| > 1): {len(unstable_eigenvalues)}'
    
    ax2.text(0.65, 0.95, stats_text, transform=ax2.transAxes,
            fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
    
    plt.tight_layout()
    # plt.savefig('floquet_analysis_2d_simple.pdf', dpi=300)
    plt.savefig('floquet_analysis_2d_simple_4field.png', dpi=300)
    
    print(f"\nAnalysis completed:")
    print(f"  Total eigenvalues analyzed: {len(eigenvalues)}")
    print(f"  Stable eigenvalues (|ρ| < 1): {len(stable_eigenvalues)}")
    print(f"  Stable eigenvalues (|ρ| = 1): {len(one_eigenvalues)}")
    print(f"  Unstable eigenvalues (|ρ| > 1): {len(unstable_eigenvalues)}")
    if len(unstable_eigenvalues) > 0:
        print(f"  Maximum modulus: {np.max(moduli):.3f}")
    
    plt.show()
# ================================
# 4. Main Function
# ================================
def main():
    print("=" * 60)
    print("Floquet Stability Analysis in 2D Polar Coordinates")
    print(f"Parameters: μ=0.5, λ=1.78, λ₄=-2.56, ω=0.45")
    print("=" * 60)
    
    # Initialize system
    print("\nInitializing system in 2D polar coordinates...")
    system = OscillonSystem()
    
    # Plot background fields
    print("\nPlotting background fields...")
    plot_background_fields(system)
    
    # Create analyzer and compute eigenvalues
    print("\nCreating Floquet analyzer...")
    analyzer = FloquetAnalyzer(system)
    
    print("Calculating eigenvalues...")
    # eigenvalues = analyzer.compute_stability(num_modes=100)
    results = analyzer.analyze_stability_high_precision()
    
    # Create paper-format plots
    print("\nGenerating plots...")
    create_paper_format_plots(system, results)

# ================================
# 5. Run Program
# ================================
if __name__ == "__main__":
    main()





