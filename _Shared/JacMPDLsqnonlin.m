function [jac] = JacMPDLsqnonlin(structModel, expModes, simModes, ...
    eigFreqOpt, normOpt, objOpt)
% function [jac] = JacMPDLsqnonlin(structModel, expModes, simModes, ...
%     eigFreqOpt, normOpt, objOpt) 
%
%   Yang Wang, Xinjun Dong, Dan Li
%   School of Civil and Environmental Engineering
%   Georgia Institute of Technology
%   2018
%
% Revision: 1.0
%
% For implementation with MATLAB lsqnonlin, this function calculates
% the Jacobian matrix for various forms of the modal property
% difference approaches.
%
% Input:
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
%     model updating:
%       lambdaExp (n_modes x 1) - experimental eigenvalue. Here n_modes
%          refers to the number of experimental modes available
%       psiExp (n_meas x n_modes) - experimental mode shape vector at
%          measured DOFs.Here n_meas refers to the number of measured DOFs
%       measDOFs (n_meas x 1) - measured DOFs
%       lambdaWeights (n_modes x 1) - weighting factor for eigenvalue
%       psiWeights (n_modes x 1) - weighting factor for eigenvector
%
%   simModes - a structure array with simulated modal properties for
%     model updating:
%       Lambda (n_modes x 1) - simulated eigenvalue
%       psi_m  (n_meas x n_modes) - simulated mode shape vector at
%          measured DOFs
%       psi    (N x n_modes) - simulated mode shape vector at all DOFs

%   eigFreqOpt
%       0 - eigenvalue difference
%       1 - angular frequency difference (rad/s)
%       2 - ordinary frequency differnce (Hz)
%
%   normOpt
%       1 - normalize eigenvecotr qi-th entry equal to 1
%       2 - normalize eigenvecotr norm equal to 1
%
%   objOpt
%       1 - MAC value formulation
%       2 - eigenvector difference formulation
%
% Output:
%   jac: the Jacobian of the objective function (matrix)

n_alpha = length( structModel.K_j ) ;
N = size( structModel.K0, 1 );

omegaSim = sqrt( simModes.Lambda );
omegaExp = sqrt( expModes.lambdaExp );

n_meas = expModes.n_meas;
n_modes = expModes.n_modes;

for i = 1 : n_modes
    if (normOpt == 1)
        simModes.psi_m(:,i) = simModes.psi_m(:,i) / simModes.psi_m(expModes.q(i), i);
        simModes.psi(:,i) = simModes.psi(:,i) / simModes.psi(expModes.q(i), i);
        
    elseif (normOpt == 2)
        simModes.psi_m(:,i) = simModes.psi_m(:,i) / norm( simModes.psi(:, i) );
        simModes.psi(:,i) = simModes.psi(:,i) / norm( simModes.psi(:, i) );
        
    else
        error('\nWrong nomalization option (normOpt) input for JacMPDLsqnonlin.');
    end
end

modalMass = zeros( n_modes, 1 );
for i = 1 : n_modes
    modalMass(i) = simModes.psi(:, i)' * structModel.M0 * simModes.psi(:, i);
end

% dLambda will store values for d_lambda_i / d_alpha_j
dLambda = zeros( n_modes, n_alpha );

% d_ri_eigFreqTerm will store the first part of d_r_i / d_alpha, i.e.
% the part involving the derivative of eigenvalue | angular
% frequency " ordinary frequency over alpha. For example, when eiegenvalue
% is used, the formulation is:
%       -weight_Lambda_i * D_alpha(Lambda_i) / Lambda_i^EXP
d_ri_eigFreqTerm = zeros( n_modes, n_alpha );

for i = 1 : n_modes
    for j = 1 : n_alpha
        dLambda(i,j) = simModes.psi(:,i)' * structModel.K_j{j} *...
            simModes.psi(:,i) / modalMass(i);
        if eigFreqOpt == 0
            d_ri_eigFreqTerm(i,j) = - dLambda(i,j) * expModes.lambdaWeights(i) ...
                / expModes.lambdaExp(i) ;
        elseif eigFreqOpt == 1 || eigFreqOpt == 2
            d_ri_eigFreqTerm(i,j) = - dLambda(i,j) * expModes.lambdaWeights(i) ...
                / (omegaExp(i) * 2 * omegaSim(i));
        end
    end
end

% The 3rd dimension corresponds to i in Psi_i^m -- the mode index
dPsi_m = zeros(n_meas, n_alpha, n_modes);
for i = 1 : n_modes
    for j = 1 : n_alpha
        dPsi_dAlpha_j = zeros(N, 1);
        B = structModel.K - simModes.Lambda(i) * structModel.M0;
        b = dLambda(i,j) * structModel.M0 * simModes.psi(:, i) -...
            structModel.K_j{j} * simModes.psi(:, i);
        
        if (normOpt == 1)
            % The maximum entry of Psi_m is normalized to 1
            P_i = setdiff(1 : N, expModes.q(i));
            % Cross out q_i-th row in B and b, and q_i-th column in B
            B = B(P_i, P_i);
            b = b(P_i);
            dPsi_dAlpha_j(P_i) = B \ b;
            Q_i = P_i(1 : n_meas - 1);
            % The q_i-th entry of dPsi_m remains as 0 from initialization
            dPsi_m(Q_i, j, i) = dPsi_dAlpha_j(Q_i);
            
        else
            % Length of Psi is normalized to 1.
            [~, index] = max(abs( simModes.psi(:,i) ));
            B(:, index) = zeros(N, 1);
            B(index, :) = zeros(1, N);
            B(index, index) = 1;
            b(index) = 0;
            v = B \ b;
            % R. B. Nelson, "Simplified calculation of eigenvector
            % derivatives," AIAA journal, vol. 14, pp. 1201-1205, 1976.
            c = -simModes.psi(:,i)' *  v ;
            dPsi_dAlpha_j = v + c * simModes.psi(:,i);
            dPsi_m(:, j, i) = dPsi_dAlpha_j(1 : n_meas);
        end
    end
end

if (objOpt == 1) 
    jac = zeros(2 * n_modes, n_alpha);
else
    if(normOpt == 1)
        jac = zeros(n_meas * n_modes, n_alpha);
    else
        jac = zeros((n_meas + 1) * n_modes, n_alpha);
    end
end

for i = 1 : n_modes
    psiExp = expModes.psiExp(:,i);
    psiSim = simModes.psi_m(:,i);
    
    if (objOpt == 1)
        % For the MAC value formulation, d_ri_MAC is the second part of
        % d_r_i/d_alpha that involves MAC value:
        MACValue = MAC(psiExp, psiSim);
        d_ri_MAC = (-expModes.psiWeights(i) / sqrt(MACValue)) * ...
            (psiExp' / (psiExp' * psiSim) - psiSim' / norm(psiSim)^2) ...
            * dPsi_m(:, :, i);
        jac((i - 1) * 2 + 1 : i * 2, :) = [d_ri_eigFreqTerm(i,:); d_ri_MAC];
    else
        % For the eigenvector difference formulation, d_ri_psi is the
        % second part of d_r_i/d_alpha that involves eigenvectors:
        %      -Qi * D_alpha(Psi_i^m) * weight_Psi_i
        if(normOpt == 1)
            Q_i = setdiff(1 : n_meas, expModes.q(i));
            d_ri_psi = -dPsi_m(Q_i, :, i) * expModes.psiWeights(i);
            jac((i - 1) * n_meas + 1 : i * n_meas, :) = ...
                [d_ri_eigFreqTerm(i,:); d_ri_psi];
        else
            d_ri_psi = -expModes.psiWeights(i) * ...
                ((psiExp' * psiSim* eye(n_meas) + psiSim *...
                psiExp') / norm(psiSim)^2 - 2 * (psiSim * psiExp' * (psiSim * psiSim'))...
                / norm(psiSim)^4) * dPsi_m(:, :, i);
            jac((i - 1) * (n_meas + 1) + 1 : i * (n_meas + 1),:) = ...
                [d_ri_eigFreqTerm(i,:); d_ri_psi];
        end
    end
end
