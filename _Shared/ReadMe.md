# ReadMe.md

**Table of Contents**

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

# StructModelUpdating.m
This function performs one run of the finite element model updating usingmfrequency domain modal properties. The starting point for this run can be provided as optimzOpts.x0.
## Syntax
function updtResults = StructModelUpdating (structModel, expModes)
function updtResults = StructModelUpdating (structModel, expModes, updatingOpts)
function updtResults = StructModelUpdating (structModel, expModes, [], optimzOpts)
function updtResults = StructModelUpdating (structModel, expModes, updatingOpts, optimzOpts)
## Description
### Input Arguments
#### structModel - a MATLAB structure array with following fields of structural model information:
|Field Name    |Dimension        |Description                    |
| ------------ | --------------- | ------------------------------|
|M0            |N x N            |mass matrix (assumed accurate enough and no need to update in current revision). Here N refers to the number of degrees of freedom of the finite element model|
|K0            |N x N            |nominal stiffness matrix constructed with nominal parameter values|
|K_j           |N x N x n_alpha  |influence matrix corresponding to updating variables (Note: the third dimension of K_j should be equal to the number of updating variables). Here n_alpha refers the number of stiffness updating variables|
#### expModes - a MATLAB structure array with experimental modal properties for model updating:
|Field Name    |Dimension        |Description                    |
| ------------ | --------------- | ------------------------------|
|lambdaExp     |n_modes x 1      |experimental eigenvalue. Here n_modes refers to the number of experimental modes available|
|psiExp        |n_meas x n_modes |experimental mode shape vector at measured DOFs. Here n_meas refers to the number of measured DOFs|
|measDOFs      |n_meas x 1       |measured DOFs|
|lambdaWeights |n_modes x 1      |weighting factor for eigenvalue|
|psiWeights    |n_modes x 1      |weighting factor for eigenvector|
#### updatingOpts - a MATLAB structure array with model updating options:
|Field Name    |Description                    |
| ------------ | ------------------------------|
|formID        |formulation ID number (default: 1)
- 1: Case 1 - conventional modal property difference formulation using MAC values
- 2: Case 2 - modal property difference formulation with eigenvector difference formulation|
|modeMatch     |Option for the matching method between simulated and experimental modes (default: 1)
- 1: Match by the MAC value between the pair of simulated and experimental mode shape vectors.
- 2: Strictly match the first designated simulated mode with the first experimental mode (see simModesForExpMatch on how to designate the simulated modes), the second designated simulated mode with the second experimental mode, etc. Only use this option when we are confident all the lowest few experimental modes are captured, i.e. there is no missing or unmeasured/undetected mode from the experimental data.|
|simModesForExpMatch|designate simulated modes obtained from FE model for matching with experimental modes
- If modeMatch = 1, Set simModesForExpMatch as an integer representing the number of simulated modes that will be compared with experimental modes for similarity matching by MAC value. The matched pair will be used for evaluating objective function value. (default: min(n_modes x 2, N))
- If modeMatch = 2, Set simModesForExpMatch as a (n_modes x 1) array.  For evaluating the objective function, the first experimental mode will be matched with simModesForExpMatch(1)-th simulated mode; the second experimental mode will be matched with simModesForExpMatch(2)-th mode, etc.|
|x_ub |upper bounds of updating variables
- formID < 3 - n_alpha x 1 (default: [])
- formID = 3 - (n_alpha + n_unmeas x n_modes) x 1 (default: [])
Here n_unmeas refers to the number of unmeasured DOFs|
|x_lb |lowe bounds of updating variables
- formID < 3 - n_alpha x 1 (default: [])
- formID = 3 - (n_alpha + n_unmeas x n_modes) x 1 (default: [])|

WARNING: when using Levenberg-Marquardt optimization algorithm in MATLAB, setting the upper and lower bounds of updating variables has no effect because the MATLAB L-M implementation does not accept bounds. The optimization may provide an infeasible out-of-the-bound solution. The user needs to verify the feasibility of the solution.
#### optimzOpts - optimization options. The current revision supports MATLAB lsqnonlin function.
|Field Name    |Description                    |
| ------------ | ------------------------------|
|maxIter       |maximum iterations of optimization process (default: 400)|
|maxFunEvals   |maximum number of function evaluations allowed (default: 100 x n_alpha)|
|tolFun        |termination tolerance on the value change of the objective function between two iterations (default: 1e-6)|
|tolX          |termination tolerance on the value change of the optimization vector variables between two iterations (default: 1e-6)|
|gradSel       |selection for gradient calculation
- 'on': calculate search gradient with user-defined Jacobian matrix. With r representing the residual vector whose length square is the objective function value, The gradient [d_f / d_alpha]' = [(d_f / d_r) * ((d_r / d_alpha)]'
- 'off': let MATLAB numerically calculate gradient matrix using finite difference method (default)|
|optAlgorithm  |optimization algorithm
- 'trust-region-reflective' algorithm
- 'Levenberg-Marquardt' algorithm (default) |
|x0            |initial value for updating variables
- formID < 3 - n_alpha x 1 
- formID = 3 - (n_alpha + n_unmeas x n_modes) x 1 
(default: zero vector)|
### Output Arguments
#### updtResults - A structure array with model updating results:
|Field Name    |Description                    |
| ------------ | ------------------------------|
|xOpt          |optimal value of updating variables |
|fvalOpt       |optimal value of objective function value|
|exitFlag      |exit flag of MATLAB lsqnonlin (check MATLAB lsqnolin help for detail: https://www.mathworks.com/help/optim/ug/lsqnonlin.html) |
|gradient      |gradient of objective function at alphaOpt|

# ModelUpdatingObjective.m
This function calculates the objective residuals of the optimization problem for model updating. When MATLAB lsqnonlin solver is used, the output contains a vector (r) whose entries are the residuals.
## Syntax
r = ModelUpdatingObjective(alpha, structModel, expModes, simModes, updatingOpts)
## Description
### Input Arguments
#### x - a vector with values of the optimization variables
#### structModel - a MATLAB structure array with following fields of structural model information:
|Field Name    |Dimension        |Description                    |
| ------------ | --------------- | ------------------------------|
|M0            |N x N            |mass matrix (assumed accurate enough and no need to update in current revision). Here N refers to the number of degrees of freedom of the finite element model|
|K0            |N x N            |nominal stiffness matrix constructed with nominal parameter values|
|K_j           |N x N x n_alpha  |influence matrix corresponding to updating variables (Note: the third dimension of K_j should be equal to the number of updating variables). Here n_alpha refers the number of stiffness updating variables|
|K             |N x N            |stiffness matrix constructed with the current alpha values, using K0 and K_j|
#### expModes - a MATLAB structure array with experimental modal properties for model updating:
|Field Name    |Dimension        |Description                    |
| ------------ | --------------- | ------------------------------|
|lambdaExp     |n_modes x 1      |experimental eigenvalue. Here n_modes refers to the number of experimental modes available|
|psiExp        |n_meas x n_modes |experimental mode shape vector at measured DOFs. Here n_meas refers to the number of measured DOFs|
|measDOFs      |n_meas x 1       |measured DOFs|
|lambdaWeights |n_modes x 1      |weighting factor for eigenvalue|
|psiWeights    |n_modes x 1      |weighting factor for eigenvector|
#### simModes - a MATLAB structure array with simulated modal properties for model updating:
|Field Name |Dimension        | Description                    |
| ----------|---------------- | ------------------------------ |
|Lambda     |n_modes x 1      |simulated eigenvalue|
|psi_m      |n_meas x n_modes |simulated mode shape vector at measured DOFs|
|psi        |N x n_modes      |simulated mode shape vector at all DOFs|
#### updatingOpts - a MATLAB structure array with model updating options:
|Field Name    |Description                    |
| ------------ | ------------------------------|
|formID        |formulation ID number (default: 1)
- 1: Case 1 - conventional modal property difference formulation using MAC values
- 2: Case 2 - modal property difference formulation with eigenvector difference formulation|
### Output Arguments
#### r - the objective residual vector r(x)

# ModelUpdatingJacobian.m
This function calculates the analytical Jacobian matrix of the optimization problem for model updating. When MATLAB lsqnonlin solver is used, the Jacobian matrix is the derivative the objective residual vector (r) over updating variables alpha: d_r / d_alpha.
## Syntax
jac = ModelUpdatingJacobian(alpha, structModel, expModes, simModes, updatingOpts)
## Description
### Input Arguments
#### x - a vector with values of the optimization variables
#### structModel - a MATLAB structure array with following fields of structural model information:
|Field Name    |Dimension        |Description                    |
| ------------ | --------------- | ------------------------------|
|M0            |N x N            |mass matrix (assumed accurate enough and no need to update in current revision). Here N refers to the number of degrees of freedom of the finite element model|
|K0            |N x N            |nominal stiffness matrix constructed with nominal parameter values|
|K_j           |N x N x n_alpha  |influence matrix corresponding to updating variables (Note: the third dimension of K_j should be equal to the number of updating variables). Here n_alpha refers the number of stiffness updating variables|
|K             |N x N            |stiffness matrix constructed with the current alpha values, using K0 and K_j|
#### expModes - a MATLAB structure array with experimental modal properties for model updating:
|Field Name    |Dimension        |Description                    |
| ------------ | --------------- | ------------------------------|
|lambdaExp     |n_modes x 1      |experimental eigenvalue. Here n_modes refers to the number of experimental modes available|
|psiExp        |n_meas x n_modes |experimental mode shape vector at measured DOFs. Here n_meas refers to the number of measured DOFs|
|measDOFs      |n_meas x 1       |measured DOFs|
|lambdaWeights |n_modes x 1      |weighting factor for eigenvalue|
|psiWeights    |n_modes x 1      |weighting factor for eigenvector|
#### simModes - a MATLAB structure array with simulated modal properties for model updating:
|Field Name |Dimension        | Description                    |
| ----------|---------------- | ------------------------------ |
|Lambda     |n_modes x 1      |simulated eigenvalue|
|psi_m      |n_meas x n_modes |simulated mode shape vector at measured DOFs|
|psi        |N x n_modes      |simulated mode shape vector at all DOFs|
#### updatingOpts - a MATLAB structure array with model updating options:
|Field Name    |Description                    |
| ------------ | ------------------------------|
|formID        |formulation ID number (default: 1)
- 1: Case 1 - conventional modal property difference formulation using MAC values
- 2: Case 2 - modal property difference formulation with eigenvector difference formulation|
### Output Arguments
#### jac - the Jacobian of the objective function (matrix)

# LsqnonlinObjJac.m
For implementation with MATLAB lsqnonlin, this function calculates the objective residual vector r(x) and returns as the first argument. When the function is called with two  output arguments:
[r, jac] = LsqnonlinObjJac(alpha, structModel, expModes, updatingOpts)
The function evaluates the user-provided analytical Jacobian matrix (d_r/d_x), and returns in the second argument.
## Syntax
r = LsqnonlinObjJac(x, structModel, expModes, updatingOpts)
[r, jac] = LsqnonlinObjJac(x, structModel, expModes, updatingOpts)
## Description
### Input Arguments
#### x - a vector with values of the optimization variables
#### structModel - a MATLAB structure array with following fields of structural model information:
|Field Name    |Dimension        |Description                    |
| ------------ | --------------- | ------------------------------|
|M0            |N x N            |mass matrix (assumed accurate enough and no need to update in current revision). Here N refers to the number of degrees of freedom of the finite element model|
|K0            |N x N            |nominal stiffness matrix constructed with nominal parameter values|
|K_j           |N x N x n_alpha  |influence matrix corresponding to updating variables (Note: the third dimension of K_j should be equal to the number of updating variables). Here n_alpha refers the number of stiffness updating variables|
#### expModes - a MATLAB structure array with experimental modal properties for model updating:
|Field Name    |Dimension        |Description                    |
| ------------ | --------------- | ------------------------------|
|lambdaExp     |n_modes x 1      |experimental eigenvalue. Here n_modes refers to the number of experimental modes available|
|psiExp        |n_meas x n_modes |experimental mode shape vector at measured DOFs. Here n_meas refers to the number of measured DOFs|
|measDOFs      |n_meas x 1       |measured DOFs|
|lambdaWeights |n_modes x 1      |weighting factor for eigenvalue|
|psiWeights    |n_modes x 1      |weighting factor for eigenvector|
#### updatingOpts - a MATLAB structure array with model updating options:
|Field Name    |Description                    |
| ------------ | ------------------------------|
|formID        |formulation ID number (default: 1)
- 1: Case 1 - conventional modal property difference formulation using MAC values
- 2: Case 2 - modal property difference formulation with eigenvector difference formulation|
|modeMatch     |Option for the matching method between simulated and experimental modes (default: 1)
- 1: Match by the MAC value between the pair of simulated and experimental mode shape vectors.
- 2: Strictly match the first designated simulated mode with the first experimental mode (see simModesForExpMatch on how to designate the simulated modes), the second designated simulated mode with the second experimental mode, etc. Only use this option when we are confident all the lowest few experimental modes are captured, i.e. there is no missing or unmeasured/undetected mode from the experimental data.|
|simModesForExpMatch|designate simulated modes obtained from FE model for matching with experimental modes
- If modeMatch = 1, Set simModesForExpMatch as an integer representing the number of simulated modes that will be compared with experimental modes for similarity matching by MAC value. The matched pair will be used for evaluating objective function value. (default: min(n_modes x 2, N))
- If modeMatch = 2, Set simModesForExpMatch as a (n_modes x 1) array.  For evaluating the objective function, the first experimental mode will be matched with simModesForExpMatch(1)-th simulated mode; the second experimental mode will be matched with simModesForExpMatch(2)-th mode, etc.|
|x_ub |upper bounds of updating variables
- formID < 3 - n_alpha x 1 (default: [])
- formID = 3 - (n_alpha + n_unmeas x n_modes) x 1 (default: [])
Here n_unmeas refers to the number of unmeasured DOFs|
|x_lb |lowe bounds of updating variables
- formID < 3 - n_alpha x 1 (default: [])
- formID = 3 - (n_alpha + n_unmeas x n_modes) x 1 (default: [])|

WARNING: when using Levenberg-Marquardt optimization algorithm in MATLAB, setting the upper and lower bounds of updating variables has no effect because the MATLAB L-M implementation does not accept bounds. The optimization may provide an infeasible out-of-the-bound solution. The user needs to verify the feasibility of the solution.
### Output Arguments
#### r - the objective residual vector r(x)
#### jac - the Jacobian of the objective function (matrix)

# ObjFuncMPDLsqnonlin.m
For implementation with MATLAB lsqnonlin, this function calculates the objective residual vector for various forms of the modal property difference approaches.
## Syntax
r =  ObjFuncMPDLsqnonlin(expModes, simModes, eigFreqOpt, normOpt, objOpt)
## Description
### Input Arguments
#### expModes - a MATLAB structure array with experimental modal properties for model updating:
|Field Name    |Dimension        | Description                    |
| ------------ | --------------- | ------------------------------ |
|lambdaExp     |n_modes x 1      |experimental eigenvalue. Here n_modes refers to the number of experimental modes available |
|psiExp        |n_meas x n_modes |experimental mode shape vector at measured DOFs. Here n_meas refers to the number of measured DOFs |
|measDOFs      |n_meas x 1       |measured DOFs|
|lambdaWeights |n_modes x 1      |weighting factor for eigenvalue|
|psiWeights    |n_modes x 1      |weighting factor for eigenvector|
#### simModes - a MATLAB structure array with simulated modal properties for model updating:
|Field Name |Dimension        | Description                    |
| ----------|---------------- | ------------------------------ |
|Lambda     |n_modes x 1      |simulated eigenvalue|
|psi_m      |n_meas x n_modes |simulated mode shape vector at measured DOFs|
|psi        |N x n_modes      |simulated mode shape vector at all DOFs|
#### eigFreqOpt - eigenfrequency options
|Values | Description         |
|-------|---------------------|
|0      |use eigenvalue difference|
|1      |angular frequency difference (rad/s)|
|2      |ordinary frequency differnce (Hz)|
#### normOpt - norm options
|Values | Description         |
|-------|---------------------|
|1      |normalize the qi-th entry of the eigenvector to 1|
|2      |normalize the length of the eigenvector to 1|
#### objOpt - objective function options
|Values | Description         |
|-------|---------------------|
|1      |MAC value formulation|
|2      |eigenvector difference formulation|
### Output Arguments
#### r - the objective residual vector r(x)

# JacMPDLsqnonlin.m
For implementation with MATLAB lsqnonlin, this function calculates the Jacobian matrix for various forms of the modal property difference approaches.
## Syntax
[jac] = JacMPDLsqnonlin(structModel, expModes, simModes, eigFreqOpt, normOpt, objOpt)
## Description
### Input Arguments
#### structModel - a MATLAB structure array with following fields of structural model information:
|Field Name    |Dimension        |Description                    |
| ------------ | --------------- | ------------------------------|
|M0            |N x N            |mass matrix (assumed accurate enough and no need to update in current revision). Here N refers to the number of degrees of freedom of the finite element model|
|K0            |N x N            |nominal stiffness matrix constructed with nominal parameter values|
|K_j           |N x N x n_alpha  |influence matrix corresponding to updating variables (Note: the third dimension of K_j should be equal to the number of updating variables). Here n_alpha refers the number of stiffness updating variables|
|K             |N x N            |stiffness matrix constructed with the current alpha values, using K0 and K_j|
#### expModes - a MATLAB structure array with experimental modal properties for model updating:
|Field Name    |Dimension        |Description                    |
| ------------ | --------------- | ------------------------------|
|lambdaExp     |n_modes x 1      |experimental eigenvalue. Here n_modes refers to the number of experimental modes available|
|psiExp        |n_meas x n_modes |experimental mode shape vector at measured DOFs. Here n_meas refers to the number of measured DOFs|
|measDOFs      |n_meas x 1       |measured DOFs|
|lambdaWeights |n_modes x 1      |weighting factor for eigenvalue|
|psiWeights    |n_modes x 1      |weighting factor for eigenvector|
#### simModes - a MATLAB structure array with simulated modal properties for model updating:
|Field Name |Dimension        | Description                    |
| ----------|---------------- | ------------------------------ |
|Lambda     |n_modes x 1      |simulated eigenvalue|
|psi_m      |n_meas x n_modes |simulated mode shape vector at measured DOFs|
|psi        |N x n_modes      |simulated mode shape vector at all DOFs|
#### eigFreqOpt - eigenfrequency options
|Values | Description         |
|-------|---------------------|
|0      |use eigenvalue difference|
|1      |angular frequency difference (rad/s)|
|2      |ordinary frequency differnce (Hz)|
#### normOpt - norm options
|Values | Description         |
|-------|---------------------|
|1      |normalize the qi-th entry of the eigenvector to 1|
|2      |normalize the length of the eigenvector to 1|
#### objOpt - objective function options
|Values | Description         |
|-------|---------------------|
|1      |MAC value formulation|
|2      |eigenvector difference formulation|
### Output Arguments
#### jac - the Jacobian of the objective function (matrix)

# MAC.m
This function calculates the modal assurance criterion between vec and each column in vecs (Reference: R. J. Allemang and D. L. Brown, "A correlation coefficient for modal vector analysis," in Proceedings of the 1st international modal analysis conference, pp. 110-116, 1982.)
## Syntax
macs = MAC(vec, vecs)
## Description
### Input Arguments
#### vec - one vector with dimension n x 1
#### vecs - a matrix with dimension n x m
### Output Arguments
#### macs - the modal assurance criterion between vec and each column in vecs
