function jac = ModelUpdatingJacobian(x, structModel, expModes, ...
    simModes, updatingOpts)
% function jac = ModelUpdatingJacobian(alpha, structModel, expModes, 
%   simModes, updatingOpts) 
%
%   Yang Wang, Xinjun Dong, Dan Li, Yu Otsuki
%   School of Civil and Environmental Engineering
%   Georgia Institute of Technology
%   2018
%
% Revision: 1.1
%
% This function calculates the analytical Jacobian matrix of the
% optimization problem for model updating. When MATLAB lsqnonlin solver is
% used, the Jacobian matrix is the derivative the objective residual vector
% (r) over updating variables alpha: d_r / d_alpha.
%
% Input:
%   x - optimization variables
%   structModel - a structure array with following fields of structural
%   model information:
%       M0 (N x N)- mass matrix (assumed accurate enough and no need to
%          update in current revision). Here N refers to the number of
%          degrees of freedom of the finite element model
%       K0 (N x N) - nominal stiffness matrix constructed with nominal
%          parameter values
%       K_j {N x N x n_alpha} - influence matrix corresponding to updating
%          variables (Note: the third dimension of K_j should be
%          equal to the number of updating variables). Here n_alpha refers
%          the number of stiffness updating variables
%       K (N x N) - stiffness matrix constructed with the current alpha
%         values, using K0 and K_j
%
%   expModes - a structure array with experimental modal properties for
%   model updating:
%       lambdaExp (n_modes x 1) - experimental eigenvalue. Here n_modes
%          refers to the number of experimental modes available
%       psiExp (n_meas x n_modes) - experimental mode shape vector at
%          measured DOFs. Here n_meas refers to the number of measured DOFs 
%       measDOFs (n_meas x 1) - measured DOFs  
%       lambdaWeights (n_modes x 1) - weighting factor for eigenvalue
%       psiWeights (n_modes x 1) - weighting factor for eigenvector
%       Lambda (n_modes x 1) - simulated eigenvalue
%       psi_m  (n_meas x n_modes) - simulated mode shape vector at 
%          measured DOFs
%       psi  (N x n_modes) - simulated mode shape vector at all DOFs
%
%   simModes - a structure array with simulated modal properties for
%     model updating:
%       Lambda (n_modes x 1) - simulated eigenvalue
%       psi_m  (n_meas x n_modes) - simulated mode shape vector at
%          measured DOFs
%       psi    (N x n_modes) - simulated mode shape vector at all DOFs
%
%   updatingOpts - a structure array with model updating options:
%       formID - formulation ID number (default: 1)
%          Case 1 - conventional modal property difference formulation using MAC values
%            1.0: eigenvalue difference, normalize eigenvecotr maximum entry equal to 1
%            1.1: angular frequency difference (rad/s), normalize eigenvecotr qi-th entry equal to 1
%            1.2: ordinary frequency differnce (Hz), normalize eigenvecotr qi-th entry equal to 1
%          Case 2 - modal property difference formulation with eigenvector difference formulation
%            2.0: eigenvalue difference, normalize eigenvecotr maximum entry equal to 1
%            2.1: angular frequency difference (rad/s), normalize eigenvecotr qi-th entry equal to 1
%            2.2: ordinary frequency differnce (Hz), normalize eigenvecotr qi-th entry equal to 1
%            2.3: eigenvalue difference, normalize eigenvecotr norm equal to 1
%            2.4: angular frequency difference (rad/s), normalize eigenvecotr norm equal to 1
%            2.5: ordinary frequency differnce (Hz), normalize eigenvecotr norm equal to 1
%          Case 3 - modal dynamic residual formulation
%            3.0: eigenvalue difference, normalize eigenvecotr maximum entry equal to 1
%           
% Output:
%   jac: the Jacobian of the objective function

switch updatingOpts.formID
    case 1.0
        jac = JacMPDLsqnonlin(structModel, expModes, simModes, 0, 1, 1);
    case 1.1
        jac = JacMPDLsqnonlin(structModel, expModes, simModes, 1, 1, 1);
    case 1.2
        jac = JacMPDLsqnonlin(structModel, expModes, simModes, 2, 1, 1);
    case 2.0
        jac = JacMPDLsqnonlin(structModel, expModes, simModes, 0, 1, 2);
    case 2.1
        jac = JacMPDLsqnonlin(structModel, expModes, simModes, 1, 1, 2);
    case 2.2
        jac = JacMPDLsqnonlin(structModel, expModes, simModes, 2, 1, 2);
    case 2.3
        jac = JacMPDLsqnonlin(structModel, expModes, simModes, 0, 2, 2);
    case 2.4
        jac = JacMPDLsqnonlin(structModel, expModes, simModes, 1, 2, 2);
    case 2.5
        jac = JacMPDLsqnonlin(structModel, expModes, simModes, 2, 2, 2);
    case 3.0
        jac = JacMDRLsqnonlin(x, structModel, expModes);
    otherwise
        error('Wrong objective option!');
end
