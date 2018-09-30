# ReadMe.md

**Table of Contents**

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
|Field Name |Dimension       |Description                    |
|:----------|:---------------|:------------------------------|
|M0         |N x N           |mass matrix (assumed accurate enough and no need to update in current revision). Here N refers to the number of degrees of freedom of the finite element model |
|K0         |N x N           |nominal stiffness matrix constructed with nominal parameter values |
|K_j        |N x N x n_alpha |influence matrix corresponding to updating variables (Note: the third dimension of K_j should be equal to the number of updating variables). Here n_alpha refers the number of stiffness updating variables |
#### expModes - a MATLAB structure array with experimental modal properties for model updating:
|Field Name    |Dimension        |Description                    |
|:-------------|:----------------|:------------------------------|
|lambdaExp     |n_modes x 1      |experimental eigenvalue. Here n_modes refers to the number of experimental modes available |
|psiExp        |n_meas x n_modes |experimental mode shape vector at measured DOFs. Here n_meas refers to the number of measured DOFs |
|measDOFs      |n_meas x 1       |measured DOFs |
|lambdaWeights |n_modes x 1      |weighting factor for eigenvalue |
|psiWeights    |n_modes x 1      |weighting factor for eigenvector |
#### updatingOpts - a MATLAB structure array with model updating options:
<table>
  <tr>
    <th text-align="left">Field Name</th>
    <th text-align="left">Description</th>
  </tr>
  <tr>
    <td>formID</td>
    <td>
    	formulation ID number (default: 1) <br>
        - 1: Case 1 - conventional modal property difference formulation using MAC values <br>
        - 2: Case 2 - modal property difference formulation with eigenvector difference formulation <br>
    </td>
  </tr>
  <tr>
    <td>modeMatch</td>
    <td>
    	Option for the matching method between simulated and experimental modes (default: 1) <br>
        - 1: Match by the MAC value between the pair of simulated and experimental mode shape vectors. <br>
        - 2: Strictly match the first designated simulated mode with the first experimental mode (see simModesForExpMatch on how to designate the simulated modes), the second designated simulated mode with the second experimental mode, etc. Only use this option when we are confident all the lowest few experimental modes are captured, i.e. there is no missing or unmeasured/undetected mode from the experimental data.
    </td>
  </tr>
  <tr>
  	<td>simModesForExpMatch</td>
    <td>
    	designate simulated modes obtained from FE model for matching with experimental modes <br>
        - If modeMatch = 1, Set simModesForExpMatch as an integer representing the number of simulated modes that will be compared with experimental modes for similarity matching by MAC value. The matched pair will be used for evaluating objective function value. (default: min(n_modes x 2, N)) <br>
        - If modeMatch = 2, Set simModesForExpMatch as a (n_modes x 1) array.  For evaluating the objective function, the first experimental mode will be matched with simModesForExpMatch(1)-th simulated mode; the second experimental mode will be matched with simModesForExpMatch(2)-th mode, etc. <br>
    </td>
  </tr>
  <tr>
  	<td>x_ub</td>
    <td>
    	upper bounds of updating variables <br>
        - formID < 3 : n_alpha x 1 (default: []) <br>
        - formID = 3 : (n_alpha + n_unmeas x n_modes) x 1 (default: []) <br>
        Here n_unmeas refers to the number of unmeasured DOFs <br>
    </td>
  </tr>
  <tr>
  	<td>x_ub</td>
    <td>
    	lower bounds of updating variables <br>
        - formID < 3 : n_alpha x 1 (default: []) <br>
        - formID = 3 : (n_alpha + n_unmeas x n_modes) x 1 (default: []) <br>
    </td>
  </tr>
</table>

WARNING: when using Levenberg-Marquardt optimization algorithm in MATLAB, setting the upper and lower bounds of updating variables has no effect because the MATLAB L-M implementation does not accept bounds. The optimization may provide an infeasible out-of-the-bound solution. The user needs to verify the feasibility of the solution.
#### optimzOpts - optimization options. The current revision supports MATLAB lsqnonlin function.
<table>
  <tr>
    <th text-align="left">Field Name</th>
    <th text-align="left">Description</th>
  </tr>
  <tr>
    <td>maxIter</td>
    <td>maximum iterations of optimization process (default: 400)</td>
  </tr>
  <tr>
    <td>maxFunEvals</td>
    <td>maximum number of function evaluations allowed (default: 100 x n_alpha)</td>
  </tr>
  <tr>
  	<td>tolFun</td>
    <td>termination tolerance on the value change of the objective function between two iterations (default: 1e-6)</td>
  </tr>
  <tr>
    <td>tolX</td>
    <td>termination tolerance on the value change of the optimization vector variables between two iterations (default: 1e-6)</td>
  </tr>
  <tr>
    <td>gradSel</td>
    <td>
      selection for gradient calculation <br>
      - 'on': calculate search gradient with user-defined Jacobian matrix. With r representing the residual vector whose length square is the objective function value, The gradient [d_f / d_alpha]' = [(d_f / d_r) * ((d_r / d_alpha)]' <br>
      - 'off': let MATLAB numerically calculate gradient matrix using finite difference method (default)
    </td>
  </tr>
  <tr>
  	<td>optAlgorithm</td>
    <td>
    	optimization algorithm <br>
        - 'trust-region-reflective' algorithm <br>
        - 'Levenberg-Marquardt' algorithm (default) <br>
    </td>
  </tr>
  <tr>
  	<td>x0</td>
    <td>
    	initial value for updating variables <br>
        - formID < 3 - n_alpha x 1 <br>
        - formID = 3 - (n_alpha + n_unmeas x n_modes) x 1 <br>
        (default: zero vector) <br>
    </td>
  </tr>
</table>

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
<a href="https://www.codecogs.com/eqnedit.php?latex=\min_{\mathbf{\alpha}}&space;\sum_{i=1}^{n_\mathrm{modes}}\left&space;\{&space;(\frac{{\lambda}_i^\mathrm{EXP}-\lambda_i(\mathbf{\alpha})}{\lambda_i^\mathrm{EXP}}\cdot&space;w_{\lambda_i})^2&space;&plus;&space;(\frac{1-\sqrt{\mathrm{MAC}_i(\mathbf{\alpha})}}{\sqrt{\mathrm{MAC}_i(\mathbf{\alpha})}}\cdot&space;w_{\mathbf{\psi&space;}_i})^2&space;\right&space;\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\min_{\mathbf{\alpha}}&space;\sum_{i=1}^{n_\mathrm{modes}}\left&space;\{&space;(\frac{{\lambda}_i^\mathrm{EXP}-\lambda_i(\mathbf{\alpha})}{\lambda_i^\mathrm{EXP}}\cdot&space;w_{\lambda_i})^2&space;&plus;&space;(\frac{1-\sqrt{\mathrm{MAC}_i(\mathbf{\alpha})}}{\sqrt{\mathrm{MAC}_i(\mathbf{\alpha})}}\cdot&space;w_{\mathbf{\psi&space;}_i})^2&space;\right&space;\}" title="\min_{\mathbf{\alpha}} \sum_{i=1}^{n_\mathrm{modes}}\left \{ (\frac{{\lambda}_i^\mathrm{EXP}-\lambda_i(\mathbf{\alpha})}{\lambda_i^\mathrm{EXP}}\cdot w_{\lambda_i})^2 + (\frac{1-\sqrt{\mathrm{MAC}_i(\mathbf{\alpha})}}{\sqrt{\mathrm{MAC}_i(\mathbf{\alpha})}}\cdot w_{\mathbf{\psi }_i})^2 \right \}" /></a>
- 2: Case 2 - modal property difference formulation with eigenvector difference formulation|
<a href="https://www.codecogs.com/eqnedit.php?latex=\min_{\mathbf{\alpha}}&space;\sum_{i=1}^{n_\mathrm{modes}}\left&space;\{&space;(\frac{{\lambda}_i^\mathrm{EXP}-\lambda_i(\mathbf{\alpha})}{\lambda_i^\mathrm{EXP}}\cdot&space;w_{\lambda_i})^2&space;&plus;&space;\left&space;\|&space;\mathbf{Q}_i\left&space;\{&space;\mathbf{\psi}_i^\mathrm{EXP,m}&space;-&space;\mathbf{\psi}_i^\mathrm{m}(\mathbf{\alpha})&space;\right&space;\}\cdot&space;w_{\mathbf{\psi&space;}_i}&space;\right&space;\|_2^2&space;\right&space;\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\min_{\mathbf{\alpha}}&space;\sum_{i=1}^{n_\mathrm{modes}}\left&space;\{&space;(\frac{{\lambda}_i^\mathrm{EXP}-\lambda_i(\mathbf{\alpha})}{\lambda_i^\mathrm{EXP}}\cdot&space;w_{\lambda_i})^2&space;&plus;&space;\left&space;\|&space;\mathbf{Q}_i\left&space;\{&space;\mathbf{\psi}_i^\mathrm{EXP,m}&space;-&space;\mathbf{\psi}_i^\mathrm{m}(\mathbf{\alpha})&space;\right&space;\}\cdot&space;w_{\mathbf{\psi&space;}_i}&space;\right&space;\|_2^2&space;\right&space;\}" title="\min_{\mathbf{\alpha}} \sum_{i=1}^{n_\mathrm{modes}}\left \{ (\frac{{\lambda}_i^\mathrm{EXP}-\lambda_i(\mathbf{\alpha})}{\lambda_i^\mathrm{EXP}}\cdot w_{\lambda_i})^2 + \left \| \mathbf{Q}_i\left \{ \mathbf{\psi}_i^\mathrm{EXP,m} - \mathbf{\psi}_i^\mathrm{m}(\mathbf{\alpha}) \right \}\cdot w_{\mathbf{\psi }_i} \right \|_2^2 \right \}" /></a>
### Output Arguments
#### r - the objective residual vector r(x)

