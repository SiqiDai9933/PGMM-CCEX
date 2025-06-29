# PGMM-CCEX
R code implementation for the Penalized Generalized Method of Moments with Common Correlated Effects (PGMM-CCEX) estimator, designed for spatial panel data models with multiple structural breaks and a multifactor error structure.

The implementation consists of four R functions:

1.PGMM_CCEX.R: The main function that sources the function “initial_estimation.R”,“PGMM_ADMM.R” and “post_estimation.R” . Iterates over a grid of tuning parameter (lambda) values, identifying the optimal lambda0 that minimizes the information criterion (IC). Outputs the corresponding post-PGMM-CCEX estimates, including structural break dates, parameter estimates, t statistics, and p value.

2.initial_estimation.R: Computes the initial non-penalized GMM-CCEX estimates for each time points and generates the adaptive weights which are used in the AGFL (Adaptive Group Fused Lasso).

3.PGMM_ADMM.R: Implements the ADFL via the ADMM algorithm to detect structural break dates and their number.

4.post_estimation.R: Estimates post-PGMM-CCEX results for each period after structural breaks are identified.

Reference: Shrinkage Estimation of Spatial Panel Data Models with Multiple Structural Breaks and a Multifactor Error Structure. 2025.03
