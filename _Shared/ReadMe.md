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
<style>
th, td {
    text-align: left;
}
</style>
<table>
  <tr>
    <th>Field Name</th>
    <th>Description</th>
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
