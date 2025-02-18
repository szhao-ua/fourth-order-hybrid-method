# fourth-order-hybrid-method




## Overview

A new high-order hybrid method integrating neural networks and corrected finite differences  is developed for solving elliptic equations with irregular interfaces and discontinuous solutions. Standard fourth-order finite difference discretization becomes invalid near such interfaces due to the discontinuities and requires corrections based on Cartesian derivative jumps. In traditional numerical methods, such as the augmented matched interface and boundary (AMIB) method, these derivative jumps can be reconstructed via additional approximations and are solved together with the unknown solution in an iterative procedure. Nontrivial developments have been carried out in the AMIB method in treating sharply curved interfaces, which, however, may not work for interfaces with geometric singularities. In this work, machine learning techniques are utilized to directly predict these Cartesian derivative jumps without involving the unknown solution. To this end, physics-informed neural networks (PINNs) are trained to satisfy the jump conditions for both closed and open interfaces with possible geometric singularities. The predicted Cartesian derivative jumps can then be integrated in the corrected finite differences. The resulted discrete Laplacian can be efficiently solved by fast Poisson solvers, such as fast Fourier transform (FFT) and geometric multigrid methods, over a rectangular domain with Dirichlet boundary conditions. This hybrid method is both easy to implement and efficient. Numerical experiments in two and three dimensions demonstrate that the method achieves fourth-order accuracy for the solution and its derivatives.

<!--
This project focuses on solving the Poisson interface problem described by the following partial differential equation (PDE):
$$
\Delta u(\mathbf{x}) + \lambda(\mathbf{x}) u(\mathbf{x}) = f(\mathbf{x}), \quad \mathbf{x} \in \Omega^- \cup \Omega^+,
$$

subject to the Dirichlet boundary condition:

$$
u(\mathbf{x}) = u_b(\mathbf{x}), \quad \mathbf{x} \in \partial \Omega.
$$

The computational domain $$\ \Omega \$$ is assumed to be rectangular and is partitioned into subdomains:

$$
\Omega = \Omega^{+} \cup \Omega^{-},
$$

separated by an interface $$\ \Gamma \$$. 
The solution in each subdomain is denoted by $$\ u^{+} \$$ and $$\ u^{-} \$$, while the source terms are $$\ f^{+} \$$ and $$\ f^{-} \$$.

Across the interface $$\ \Gamma \$$, the solution exhibits jump discontinuities governed by the following conditions:

$$
[\![u(\mathbf{x})]\!] = \gamma(\mathbf{x}), \quad [\![u_n(\mathbf{x})]\!] = \rho(\mathbf{x}), \quad \mathbf{x} \in \Gamma,
$$

where $$\( [\![u(\mathbf{x})]\!] \) $$ and $$\( [\![u_n(\mathbf{x})]\!] \)$$ represent the jump in the function and its normal derivative, respectively.-->

## Features

1. **Hybrid Neural Network and Finite Difference Method**:
   - A neural network is trained to approximate the jump conditions across the interface.
   - The finite difference method enforces the governing equation in each subdomain.
   - Achieves **fourth-order accuracy** for the solution and its derivatives.



2. **Fortran Integration**:
   - Utilizes Fortran subroutines (`fft2d.f90` and `fft.f`) for efficient numerical computations.
   - The Fortran code is imported into Python as a module (`mymodule`).

---

## Requirements

The following dependencies are required to run the project:
- Python 3.8+
- PyTorch
- Functorch
- f2py (for importing Fortran subroutines)

---

## File Structure

- **`fft2d.f90`** and **`fft.f`**: Fortran subroutines for numerical computations.
- **`mymodule`**: Python module wrapping the Fortran subroutines.

---

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/szhao-ua/fourth-order-hybrid-method.git
   cd fourth-order-hybrid-method

2. Compile the Fortran subroutines using f2py:
   ```bash
   f2py -c fft2d.f90 fft.f -m mymodule


3. Run the main python script

