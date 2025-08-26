===========================================================
             GEMANOVA TOOLBOX v1.0 (MATLAB)
===========================================================

The GEMANOVA Toolbox provides functions for fitting, 
evaluating, interpreting, and validating Generalized 
Multivariate ANOVA (GEMANOVA) models in MATLAB.
Useful in datasets that follow an experimental design 
and have high-order interactions present.

-----------------------------------------------------------
FEATURES
-----------------------------------------------------------
- Tools for preprocessing (replicate averaging, cube reshaping)
- Model fitting with GEMANOVA decomposition
- Permutation testing (full-factor and pairwise)
- Variance partitioning and explained variance analysis
- Workflow scripts for fitting, evaluation, and validation

-----------------------------------------------------------
FILES
-----------------------------------------------------------
GEMANOVAprocess.m           - Pipeline for running GEMANOVA analysis
average_replicates.m        - Averages replicate measurements
MatrixToHypercube.m         - Reshapes data matrix into hypercube format
FitModel.m                  - Wrapper for fitting GEMANOVA models
VarianceExplainedPerModel.m - Computes variance explained by factor combinations
FullFactorPermutation.m     - Performs full-factor permutation tests
PairwiseFactorPermutation.m - Pairwise permutation testing of factors
gemanova.m                  - Core GEMANOVA fitting function
gemanovaInit.m              - GEMANOVA with initial estimates
pffitalpha.m                - Parameter fitting utility (alpha estimation)

-----------------------------------------------------------
INSTALLATION
-----------------------------------------------------------
1. Copy the toolbox folder to your computer.
2. Add the toolbox to the MATLAB path:
   >> addpath(genpath('path/to/GEMANOVA_toolbox'))
   >> savepath

3. Test installation by running:
   >> help gemanova

-----------------------------------------------------------
CITATION
-----------------------------------------------------------
Coded by Jokin Ezenarro, 2025. (and ChatGPT)
Based on the work of Rasmus Bro, 2002.

If you use this toolbox in your research, please cite:
- 

-----------------------------------------------------------
LICENSE
-----------------------------------------------------------
No commercial use of this software is allowed.
Otherwise, you are free to use, modify, and distribute with attribution.