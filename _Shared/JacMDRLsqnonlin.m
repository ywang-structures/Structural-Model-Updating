function [jac] = JacMDRLsqnonlin(x,structModel,expModes)
% function [jac] = JacMDRLsqnonlin(x,structModel,expModes)
%
%   Yang Wang, Xinjun Dong, Dan Li
%   School of Civil and Environmental Engineering
%   Georgia Institute of Technology
%   2018
%
% Revision: 1.1
%
% This function calculates the Jacobian matrix of the objective residual 
% vector r for the objective function of the modal dynamic residual 
% approach with respect to the optimization variables
%
% Note: the optimization variables include the stiffness parameters and 
% the unmeasured entries of the eigenvectors
%
% Input:
%   x - optimization variables, including stiffness parameters and 
%       unmeasured entries of the eigenvectors
%
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
%
%   expModes - a structure array with experimental modal properties for
%   model updating:
%       lambdaExp (n_modes x 1) - experimental eigenvalue. Here n_modes
%          refers to the number of experimental modes available
%       psiExp (n_meas x n_modes) - experimental mode shape vector at
%          measured DOFs.Here n_meas refers to the number of measured DOFs 
%       measDOFs (n_meas x 1) - measured DOFs  
%       resWeights (n_modes x 1) - weighting factor for residual
% Output:
%   jac: Jacobian matrix of the objective function

weight = expModes.resWeights;
N = size(structModel.M0,1);
n_alpha = length(structModel.K_j);
n_u = N - expModes.n_meas;

Dr_alpha = zeros(N * expModes.n_modes, n_alpha);

for i = 1:n_alpha
    for j = 1:expModes.n_modes
        psiMix(1 : expModes.n_meas,1) = expModes.psiExp(:,j);
        psiMix(expModes.n_meas + 1 : N,1) = x(n_alpha + (j-1) * n_u + 1: n_alpha + j * n_u);
        Dr_alpha((j-1) * N + 1 : j * N, i) = structModel.K_j{i} * psiMix * weight(j);
    end
end

Dr_psiU = zeros(N * expModes.n_modes, n_u * expModes.n_modes);
for i = 1 : expModes.n_modes
    K_u = structModel.K(:, expModes.n_meas + 1: end);
    M_u = structModel.M0(:, expModes.n_meas + 1: end);
    Dr_psiU((i-1) * N + 1 : i * N , (i-1) * n_u + 1: i * n_u) = (K_u - expModes.lambdaExp(i) * M_u) * weight(i);
end

jac = [Dr_alpha Dr_psiU];
