# 2QG_IM

<!-- Dynamic Badges (Update URLs) -->
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![GitHub Release](https://img.shields.io/github/v/release/yourusername/repo)](https://github.com/yourusername/repo/releases)
[![Build Status](https://img.shields.io/github/actions/workflow/status/yourusername/repo/build.yml)](https://github.com/yourusername/repo/actions)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen)](https://yourusername.github.io/repo/)

---

## Project Overview 🔍
The time evolution system governing  the fields involved in this equilibrium state constitutes the subject of quasi-geographic (QG) theory, yielding the so-called QG models which  describe  the essential characteristics of the geophysical fluid flows.  

we proposed a linearized iterative scheme(IM) to compute the  non-dimensional potential vorticity  and stream function in 2-dimensional space successively.

In practice, the Rossby number $Ro$ is small $(O(10^{-3}))$ in QG model, so a highly refined mesh is required to direct numerical simulation of model, which leads to a time-consuming calculation scheme and requires a significant amount of computational memory. To overcome this difficulty, many researchers have introduced the large eddy simulation (LES) technique.

In our work, we set the Rossby number $Ro$ not too small($Ro$ =1.5,0.088) to avoid the requirement of highly refined mesh or LES technology. Our objective is to develop an efficient algorithm that is different from traditional direct numerical simulation(DNS).


## Quick Start 🚀
### Prerequisites
- MATALB R2022a

## Run the code
- The file "Example1" and "Example2" are contain the code of Example 1 and 2,respectly in section 4 of the manuscript.
- By run the main code mentioned in the file "A_explanation", you can obtain the same results in the section 4.
- the file "Myfunction5" contains the setting parameters($Ro$, $Re$, $Fr$, $\delta$, $\sigma$) and functions($\psi$, $q$).
