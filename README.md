# R code for "Hierarchical False Discovery Rate Control for High-dimensional Survival Analysis with Interactions"

This repository contains the R implementation of the Hierarchical False Discovery Rate (HFDR) method presented in the paper titled "Hierarchical False Discovery Rate Control for High-dimensional Survival Analysis with Interactions" by Liang W., Zhang, Q., and Ma, S. (2023).

## Overview

The HFDR method is designed for high-dimensional survival analysis with interactions. It provides a flexible framework for controlling the False Discovery Rate (FDR) while accounting for the hierarchical structure of the hypotheses. This R code allows users to apply the HFDR method to their survival analysis data and effectively identify significant main effects and interactions.

## Contents

HFDR.R: This file contains the main R functions for implementing the HFDR method. Users can use these functions to perform high-dimensional survival analysis with interaction terms and control the FDR.

lasso_inference.R: In addition to the HFDR method, this repository also includes implementing the "debiased lasso" method, sourced from https://web.stanford.edu/~montanar/sslasso/code.html. 

## Usage

To use the HFDR or alternative methods, simply clone this repository and source the respective R script in your R environment. Detailed usage instructions and examples can be found within the script files.

## Citation

If you use the HFDR method or the code from this repository in your research, please cite the following paper:

Liang W., Zhang, Q., and Ma, S. (2023). Hierarchical False Discovery Rate Control for High-dimensional Survival Analysis with Interactions.

## Contact

If you have any questions, suggestions, or encounter any issues with the code, please feel free to open an issue on this GitHub repository or contact the authors.

