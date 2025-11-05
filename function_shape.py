import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# 1D 情况
print("=== 1D Case ===")
# Parameters
mu = 0.5
lamb = 1.78
lamb4 = -2.56
omega = 0.6

# Domain - 1D coordinates, x from -L to L
L = 50.0  # Domain half-length
N = 500
dx = 2 * L / N
x = np.linspace(-L, L, N+1)

# Number of variables: 2*(N+1)
def equations_1d(y):
    f0 = y[0:N+1]
    f1 = y[N+1:2*(N+1)]
    res0 = np.zeros(N+1)
    res1 = np.zeros(N+1)
    
    # For f0 - 1D Laplacian: ∇²f = f''
    for i in range(0, N+1):
        if i == 0:  # Left boundary: df/dx = 0
            laplacian_f0 = (2 * f0[1] - 2 * f0[0]) / dx**2
        elif i == N:  # Right boundary: df/dx = 0
            laplacian_f0 = (2 * f0[N-1] - 2 * f0[N]) / dx**2
        else:  # Interior points
            # 1D Laplacian: f''
            laplacian_f0 = (f0[i-1] - 2*f0[i] + f0[i+1]) / dx**2
        
        nonlinear_f0 = mu * f0[i] - 2*lamb * f0[i] * (4 * f0[i]**2 + 2 * f1[i]**2) + 2 * lamb4 * f0[i] * (-2 * f0[i]**2 + f1[i]**2)
        res0[i] = laplacian_f0 + nonlinear_f0
    
    # For f1 - special treatment to enforce both f1=0 and df1/dx=0 at boundaries
    for i in range(0, N+1):
        if i == 0:  # Left boundary: f1 = 0 and df1/dx = 0
            res1[i] = f1[0]
        elif i == 1:  # First interior point - use special treatment
            laplacian_f1 = (f1[0] - 2*f1[1] + f1[2]) / dx**2
            nonlinear_f1 = omega**2 * f1[i] + mu * f1[i] - 2*lamb * f1[i] * (4 * f0[i]**2 + 2 * f1[i]**2) - 2 * lamb4 * f1[i] * (-2 * f0[i]**2 + f1[i]**2)
            res1[i] = laplacian_f1 + nonlinear_f1
        elif i == N-1:  # Second to last point - use special treatment
            laplacian_f1 = (f1[N-2] - 2*f1[N-1] + f1[N]) / dx**2
            nonlinear_f1 = omega**2 * f1[i] + mu * f1[i] - 2*lamb * f1[i] * (4 * f0[i]**2 + 2 * f1[i]**2) - 2 * lamb4 * f1[i] * (-2 * f0[i]**2 + f1[i]**2)
            res1[i] = laplacian_f1 + nonlinear_f1
        elif i == N:  # Right boundary: f1 = 0 and df1/dx = 0
            res1[i] = f1[N]
        else:  # Interior points
            laplacian_f1 = (f1[i-1] - 2*f1[i] + f1[i+1]) / dx**2
            nonlinear_f1 = omega**2 * f1[i] + mu * f1[i] - 2*lamb * f1[i] * (4 * f0[i]**2 + 2 * f1[i]**2) - 2 * lamb4 * f1[i] * (-2 * f0[i]**2 + f1[i]**2)
            res1[i] = laplacian_f1 + nonlinear_f1
    
    return np.concatenate([res0, res1])

# Initial guess for 1D
f0_initial_1d = 0.43 * np.ones(N+1) - 0.43 * (np.cosh(x/3)**(-1))
f1_initial_1d = -0.6 * ((np.cosh(x/3))**(-1))
f1_initial_1d[0] = 0.0
f1_initial_1d[N] = 0.0
y_initial_1d = np.concatenate([f0_initial_1d, f1_initial_1d])

# Solve 1D
print("Solving 1D system...")
solution_1d = fsolve(equations_1d, y_initial_1d, maxfev=5000)
print("1D solution found!")

# Extract 1D solutions
f0_solution_1d = solution_1d[0:N+1]
f1_solution_1d = solution_1d[N+1:2*(N+1)]

# Calculate gradients for 1D
f0_grad_1d = np.gradient(f0_solution_1d, dx)
f1_grad_1d = np.gradient(f1_solution_1d, dx)

# Calculate the expression for 1D
expression_1d = (2 * (omega**2 * f1_solution_1d)**2 +
                 4 * f0_grad_1d**2 + 
                 2 * f1_grad_1d**2 -
                 mu * (2*f0_solution_1d**2 + f1_solution_1d**2) +
                 lamb * (2*f0_solution_1d**2 + f1_solution_1d**2)**2 +
                 lamb4 * (-2*f0_solution_1d**2 + f1_solution_1d**2)**2)

# 2D Polar Coordinates
print("\n=== 2D Polar Coordinates ===")
# Domain - polar coordinates, r from 0 to R
R = 50.0  # Maximum radius
N = 500
dr = R / N
r = np.linspace(0, R, N+1)

# Avoid division by zero at r=0
r_safe = r.copy()
r_safe[0] = 1e-10

def equations_2d(y):
    f0 = y[0:N+1]
    f1 = y[N+1:2*(N+1)]
    res0 = np.zeros(N+1)
    res1 = np.zeros(N+1)
    
    # For f0 - polar Laplacian: ∇²f = f'' + (1/r)f'
    for i in range(0, N+1):
        if i == 0:
            laplacian_f0 = 2 * (f0[1] - f0[0]) / dr**2
        elif i == N:
            laplacian_f0 = (2 * f0[N-1] - 2 * f0[N]) / dr**2
        else:
            second_deriv = (f0[i-1] - 2*f0[i] + f0[i+1]) / dr**2
            first_deriv = (f0[i+1] - f0[i-1]) / (2 * dr)
            laplacian_f0 = second_deriv + (1 / r_safe[i]) * first_deriv
        
        nonlinear_f0 = mu * f0[i] - 2*lamb * f0[i] * (4 * f0[i]**2 + 2 * f1[i]**2) + 2 * lamb4 * f0[i] * (-2 * f0[i]**2 + f1[i]**2)
        res0[i] = laplacian_f0 + nonlinear_f0
    
    # For f1 - polar Laplacian
    for i in range(0, N+1):
        if i == 0:
            laplacian_f1 = 2 * (f1[1] - f1[0]) / dr**2
        elif i == N:
            laplacian_f1 = (2 * f1[N-1] - 2 * f1[N]) / dr**2
        else:
            second_deriv = (f1[i-1] - 2*f1[i] + f1[i+1]) / dr**2
            first_deriv = (f1[i+1] - f1[i-1]) / (2 * dr)
            laplacian_f1 = second_deriv + (1 / r_safe[i]) * first_deriv
        
        nonlinear_f1 = omega**2 * f1[i] + mu * f1[i] - 2*lamb * f1[i] * (4 * f0[i]**2 + 2 * f1[i]**2) - 2 * lamb4 * f1[i] * (-2 * f0[i]**2 + f1[i]**2)
        res1[i] = laplacian_f1 + nonlinear_f1
    
    return np.concatenate([res0, res1])

# Initial guess for 2D
f0_initial_2d = 0.43 * np.ones(N+1) - 0.43 * (np.cosh(r/3)**(-1))
f1_initial_2d = -0.6 * ((np.cosh(r/3))**(-1))
f1_initial_2d[N] = 0.0
y_initial_2d = np.concatenate([f0_initial_2d, f1_initial_2d])

# Solve 2D
print("Solving 2D system...")
solution_2d = fsolve(equations_2d, y_initial_2d, maxfev=5000)
print("2D solution found!")

# Extract 2D solutions
f0_solution_2d = solution_2d[0:N+1]
f1_solution_2d = solution_2d[N+1:2*(N+1)]


expression_data = np.column_stack((f0_solution_2d, f1_solution_2d))
header_expression = "r expression_value"
np.savetxt('expression_2D_polar.txt', expression_data, header=header_expression, fmt='%.8f', comments='')
print("Expression data saved to 'expression_2D_polar.txt'")

# Calculate gradients for 2D
f0_grad_2d = np.gradient(f0_solution_2d, dr)
f1_grad_2d = np.gradient(f1_solution_2d, dr)

# Calculate the expression for 2D
expression_2d = (2 * (omega**2 * f1_solution_2d)**2 +
                 4 * f0_grad_2d**2 + 
                 2 * f1_grad_2d**2 -
                 mu * (2*f0_solution_2d**2 + f1_solution_2d**2) +
                 lamb * (2*f0_solution_2d**2 + f1_solution_2d**2)**2 +
                 lamb4 * (-2*f0_solution_2d**2 + f1_solution_2d**2)**2)

# 3D Spherical Coordinates
print("\n=== 3D Spherical Coordinates ===")
def equations_3d(y):
    f0 = y[0:N+1]
    f1 = y[N+1:2*(N+1)]
    res0 = np.zeros(N+1)
    res1 = np.zeros(N+1)
    
    # For f0 - spherical Laplacian: ∇²f = f'' + (2/r)f'
    for i in range(0, N+1):
        if i == 0:
            laplacian_f0 = 3 * (f0[1] - f0[0]) / dr**2
        elif i == N:
            laplacian_f0 = (2 * f0[N-1] - 2 * f0[N]) / dr**2
        else:
            second_deriv = (f0[i-1] - 2*f0[i] + f0[i+1]) / dr**2
            first_deriv = (f0[i+1] - f0[i-1]) / (2 * dr)
            laplacian_f0 = second_deriv + (2 / r_safe[i]) * first_deriv
        
        nonlinear_f0 = mu * f0[i] - 2*lamb * f0[i] * (4 * f0[i]**2 + 2 * f1[i]**2) + 2 * lamb4 * f0[i] * (-2 * f0[i]**2 + f1[i]**2)
        res0[i] = laplacian_f0 + nonlinear_f0
    
    # For f1 - spherical Laplacian
    for i in range(0, N+1):
        if i == 0:
            laplacian_f1 = 3 * (f1[1] - f1[0]) / dr**2
        elif i == N:
            laplacian_f1 = (2 * f1[N-1] - 2 * f1[N]) / dr**2
        else:
            second_deriv = (f1[i-1] - 2*f1[i] + f1[i+1]) / dr**2
            first_deriv = (f1[i+1] - f1[i-1]) / (2 * dr)
            laplacian_f1 = second_deriv + (2 / r_safe[i]) * first_deriv
        
        nonlinear_f1 = omega**2 * f1[i] + mu * f1[i] - 2*lamb * f1[i] * (4 * f0[i]**2 + 2 * f1[i]**2) - 2 * lamb4 * f1[i] * (-2 * f0[i]**2 + f1[i]**2)
        res1[i] = laplacian_f1 + nonlinear_f1
    
    return np.concatenate([res0, res1])

# Initial guess for 3D
f0_initial_3d = 0.43 * np.ones(N+1) - 0.43 * (np.cosh(r/3)**(-1))
f1_initial_3d = -0.6 * ((np.cosh(r/3))**(-1))
f1_initial_3d[N] = 0.0
y_initial_3d = np.concatenate([f0_initial_3d, f1_initial_3d])

# Solve 3D
print("Solving 3D system...")
solution_3d = fsolve(equations_3d, y_initial_3d, maxfev=5000)
print("3D solution found!")

# Extract 3D solutions
f0_solution_3d = solution_3d[0:N+1]
f1_solution_3d = solution_3d[N+1:2*(N+1)]

# Calculate gradients for 3D
f0_grad_3d = np.gradient(f0_solution_3d, dr)
f1_grad_3d = np.gradient(f1_solution_3d, dr)

# Calculate the expression for 3D
expression_3d = (2 * (omega**2 * f1_solution_3d)**2 +
                 4 * f0_grad_3d**2 + 
                 2 * f1_grad_3d**2 -
                 mu * (2*f0_solution_3d**2 + f1_solution_3d**2) +
                 lamb * (2*f0_solution_3d**2 + f1_solution_3d**2)**2 +
                 lamb4 * (-2*f0_solution_3d**2 + f1_solution_3d**2)**2)

# Plot the expressions for all three cases
plt.figure(figsize=(15, 10))

# 1D plot
plt.subplot(2, 3, 1)
plt.plot(x, f0_solution_1d, 'b-', label='f0')
plt.plot(x, f1_solution_1d, 'r-', label='f1')
plt.title('1D: f0 and f1')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.grid(True)
plt.legend()

plt.subplot(2, 3, 4)
plt.plot(x, expression_1d, 'g-')
plt.title('1D: Energy Density')
plt.xlabel('x')
plt.ylabel('Energy Density')
plt.grid(True)

# 2D plot
plt.subplot(2, 3, 2)
plt.plot(r, f0_solution_2d, 'b-', label='f0')
plt.plot(r, f1_solution_2d, 'r-', label='f1')
plt.title('2D Polar: f0 and f1')
plt.xlabel('r')
plt.ylabel('f(r)')
plt.grid(True)
plt.legend()

plt.subplot(2, 3, 5)
plt.plot(r, expression_2d, 'g-')
plt.title('2D Polar: Energy Density')
plt.xlabel('r')
plt.ylabel('Energy Density')
plt.grid(True)

# 3D plot
plt.subplot(2, 3, 3)
plt.plot(r, f0_solution_3d, 'b-', label='f0')
plt.plot(r, f1_solution_3d, 'r-', label='f1')
plt.title('3D Spherical: f0 and f1')
plt.xlabel('r')
plt.ylabel('f(r)')
plt.grid(True)
plt.legend()

plt.subplot(2, 3, 6)
plt.plot(r, expression_3d, 'g-')
plt.title('3D Spherical: Energy Density')
plt.xlabel('r')
plt.ylabel('Energy Density')
plt.grid(True)

plt.tight_layout()
plt.show()

# Create a combined plot for comparison
plt.figure(figsize=(12, 8))

plt.subplot(2, 1, 1)
plt.plot(x, expression_1d, 'b-', label='1D')
plt.plot(r, expression_2d, 'r-', label='2D Polar')
plt.plot(r, expression_3d, 'g-', label='3D Spherical')
plt.title('Comparison of Energy Density Across Dimensions')
plt.xlabel('Coordinate')
plt.ylabel('Energy Density')
plt.grid(True)
plt.legend()

plt.subplot(2, 1, 2)
# Normalize for better comparison
expression_1d_norm = expression_1d / np.max(np.abs(expression_1d))
expression_2d_norm = expression_2d / np.max(np.abs(expression_2d))
expression_3d_norm = expression_3d / np.max(np.abs(expression_3d))

plt.plot(x, expression_1d_norm, 'b-', label='1D (normalized)')
plt.plot(r, expression_2d_norm, 'r-', label='2D Polar (normalized)')
plt.plot(r, expression_3d_norm, 'g-', label='3D Spherical (normalized)')
plt.title('Normalized Energy Density Comparison')
plt.xlabel('Coordinate')
plt.ylabel('Normalized Energy Density')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

# Print some statistics
print("\n=== Expression Statistics ===")
print(f"1D - Min: {np.min(expression_1d):.6f}, Max: {np.max(expression_1d):.6f}, Mean: {np.mean(expression_1d):.6f}")
print(f"2D - Min: {np.min(expression_2d):.6f}, Max: {np.max(expression_2d):.6f}, Mean: {np.mean(expression_2d):.6f}")
print(f"3D - Min: {np.min(expression_3d):.6f}, Max: {np.max(expression_3d):.6f}, Mean: {np.mean(expression_3d):.6f}")