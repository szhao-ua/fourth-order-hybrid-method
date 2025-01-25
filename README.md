# corrected-hybrid-method




## Overview
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

separated by an embedded interface $$\ \Gamma \$$. 
<!--The solution in each subdomain is denoted by $$\ u^{+} \$$ and $$\ u^{-} \$$, while the source terms are $$\ f^{+} \$$ and $$\ f^{-} \$$.-->

Across the interface $$\ \Gamma \$$, the solution exhibits jump discontinuities governed by the following conditions:

$$
[\![u(\mathbf{x})]\!] = \gamma(\mathbf{x}), \quad [\![u_n(\mathbf{x})]\!] = \rho(\mathbf{x}), \quad \mathbf{x} \in \Gamma,
$$

where $$\( [\![u(\mathbf{x})]\!] \) $$ and $$\( [\![u_n(\mathbf{x})]\!] \)$$ represent the jump in the function and its normal derivative, respectively.

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
- NumPy
- PyTorch (for neural network training)
- f2py (for importing Fortran subroutines)

---

## File Structure

- **`fft2d.f90`** and **`fft.f`**: Fortran subroutines for numerical computations.
- **`mymodule`**: Python module wrapping the Fortran subroutines.

---

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/szhao-ua/corrected-hybrid-method.git
   cd corrected-hybrid-method

2. Compile the Fortran subroutines using f2py:
   ```bash
   f2py -c fft2d.f90 fft.f -m mymodule


3. Run the main python script

