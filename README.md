# Fractional Temporal SEIR Model with Signorini Boundary Conditions

**Authors:** Omar Elamraoui, EL Hassan Essoufi, Abderrahim Zafrar  
**© 2025 O. Elamraoui, E-H. Essoufi, A. Zafrar. All rights reserved.**

This repository contains MATLAB code for simulating, analyzing, and visualizing a **fractional-order SEIR epidemic model** with memory effects and **Signorini boundary conditions**, as presented in the paper:

> *Switching Cases for Fractional Time SEIR Model with Memory and Space Diffusion*

---

## 📁 Files Included

### 🔷 Main Simulation and Analysis Files
- **`Main_SEIR_Frac_Memory.m`**  
  Main driver script to simulate the SEIR model with fractional memory and Signorini boundary conditions.  
  Generates plots and supports parametric studies.

- **`show_Sensitivity_Param.m`**  
  Sensitivity analysis of SEIR dynamics with respect to a chosen parameter.  
  Produces time-series plots for S, E, I, R components.

- **`SEIR_Show_var_alpha.m`**  
  Studies the impact of varying the fractional-order parameter `alpha` on the SEIR dynamics.

- **`SEIR_var_alpha_fun.m`**  
  Supporting function for `SEIR_Show_var_alpha.m`.  
  Returns SEIR evolution for a set of `alpha` values.

### 🔷 Core Functions
- **`analyze_param_impact.m`**  
  Executes SEIR simulations across multiple values of parameters (`β`, `γ`, `α`, etc.).  
  Returns S, E, I, R trajectories for comparison.

- **`run_SEIR_simulation.m`**  
  Core solver using the Finite Element Method (FEM) and Uzawa’s method to simulate the SEIR model with fractional time memory and Signorini boundary conditions.

- **`UzawaSignoriniSolver.m`**  
  Regularized Uzawa solver for the variational inequality with Signorini boundary conditions.  
  Returns numerical solution, projections, multipliers, and convergence metrics.

### 🔷 Helper Functions
- **`kpde2dumsh.m`, `kpde2dstf.m`, `kpde2drhs.m`, `kpde2dmss.m`**  
  Helper functions for mesh generation, finite element assembly, and right-hand side computations.

---

## 💻 Requirements

- MATLAB **R2019b or newer** (recommended)
- No specialized toolboxes required (uses standard MATLAB functions)

---

## 🚀 How to Use

1. Run **`Main_SEIR_Frac_Memory.m`** or **`show_Sensitivity_Param.m`** to simulate the SEIR model or analyze parameter sensitivity.
2. Use **`SEIR_Show_var_alpha.m`** to explore the effect of the fractional-order parameter `alpha`.
3. Output plots will be generated and saved in `.jpg` and `.eps` formats with descriptive filenames.
4. Modify parameter ranges and simulation settings directly in the scripts for custom analysis.

---

## 📚 Acknowledgments

Finite element matrix assembly routines adapted from:

**Koko, J. (2016)** – *Fast finite element assembly in MATLAB and Octave*  
[https://perso.isima.fr/~jokoko/codes.html](https://perso.isima.fr/~jokoko/codes.html)

---

## 📬 Contact

For questions, suggestions, or further information, please contact:

**Omar Elamraoui**  
📧 [oelamraoui34@gmail.com](mailto:oelamraoui34@gmail.com)

---

