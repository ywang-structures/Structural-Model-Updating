# ReadMe.md

**Table of Contents**

[TOC]

# MultiRunModelUpdating.m
The script can perform multiple runs of model updating by starting the search from randomized starting points.
## Syntax
MultiRunModelUpdating
## Description
Before running this script, the workspace must have following variables before running the script. See documentation of StructModelUpdating.m for descriptions.

**structModel, expModes, updatingOpts, optimzOpts, filename**

If following variables do not exist, default actions will be assigned.

**numRuns**: number of starting points.  Default = 1
**randSeed**: random seed value for rng.  By default the random seed will NOT be fixed.

If multiple runs are requested through "numRuns" variable and the runs are disrupted half way, the script can pick up and continue the runs by loading previous run results from the .mat named by the "filename" variable. **IF A USER WANTS TO START ALL THE RUNS AFRESH, THE .MAT FILE MUST BE DELETED FIRST.**

# ObjFuncMPDLsqnonlin.m
For implementation with MATLAB lsqnonlin, this function calculates the objective residual vector for various forms of the modal property difference approaches.
## Syntax
r =  ObjFuncMPDLsqnonlin(expModes, simModes, eigFreqOpt, normOpt, objOpt)
## Description
### Input Arguments
#### expModes - a MATLAB structure array with experimental modal properties for model updating
|Fields           |Dimension  | Description                    |
| -------------|-------------| ------------------------------ |
| lambdaExp |$$n_modes x 1$$ | experimental eigenvalue. Here n_modes refers to the number of experimental modes available |
| psiExp         |n_meas x n_modes   |experimental mode shape vector at measured DOFs. Here n_meas refers to the number of measured DOFs     |
|measDOFs   |n_meas x 1  |measured DOFs|
|lambdaWeights |n_modes x 1 |weighting factor for eigenvalue|
|psiWeights  |n_modes x 1| weighting factor for eigenvector|


simModes - a MATLAB structure array with simulated modal properties for model updating:
Lambda (n_modes x 1) - simulated eigenvalue
psi_m  (n_meas x n_modes) - simulated mode shape vector at measured DOFs
psi    (N x n_modes) - simulated mode shape vector at all DOFs

eigFreqOpt
0 - use eigenvalue difference
1 - angular frequency difference (rad/s)
2 - ordinary frequency differnce (Hz)

normOpt
1 - normalize the qi-th entry of the eigenvector to 1
2 - normalize the length of the eigenvector to 1

objOpt
1 - MAC value formulation
2 - eigenvector difference formulation


