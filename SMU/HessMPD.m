function [hess] = HessMPD(x, structModel, expModes, updatingOpts)
% function [hess] = HesMPD(x, structModel, expModes, updatingOpts)
%
%   Yang Wang, Xinjun Dong, Dan Li, Yu Otsuki
%   School of Civil and Environmental Engineering
%   Georgia Institute of Technology
%   2018
%
% Revision: 1.0
%
% This function calculates the Hessian matrix for various forms of the 
% modal property difference approaches.
%
% Input:
%   x - a vector with values of the optimization variables
%
%   structModel - a MATLAB structure array with following fields of
%   structural model information:
%       M0 (N x N)- mass matrix (assumed accurate enough and no need to
%          update in current revision). Here N refers to the number of
%          degrees of freedom of the finite element model
%       K0 (N x N) - nominal stiffness matrix constructed with nominal
%          parameter values
%       K_j (N x N x n_alpha) - influence matrix corresponding to updating
%          variables (Note: the third dimension of K_j should be
%          equal to the number of updating variables). Here n_alpha refers
%          the number of stiffness updating variables
%       K (N x N) - stiffness matrix constructed with the current alpha
%         values, using K0 and K_j
%
%   expModes - a MATLAB structure array with experimental modal properties
%   for model updating:
%       lambdaExp (n_modes x 1) - experimental eigenvalue. Here n_modes
%          refers to the number of experimental modes available
%       psiExp (n_meas x n_modes) - experimental mode shape vector at
%          measured DOFs. Here n_meas refers to the number of measured DOFs
%       measDOFs (n_meas x 1) - measured DOFs
%       lambdaWeights (n_modes x 1) - weighting factor for eigenvalue
%       psiWeights (n_modes x 1) - weighting factor for eigenvector
%
%   updatingOpts - a MATLAB structure array with model updating options:
%       formID - formulation ID number (default: 1)
%           1: Case 1 - conventional modal property difference
%              formulation using MAC values
%           2: Case 2 - modal property difference formulation with
%             eigenvector difference formulation
%
% Output:
%   hess: the Hessian matrix of the objective function



%% Options
switch updatingOpts.formID
    case 1.0
         eigFreqOpt = 0; normOpt = 1; objOpt = 1;
    case 1.1
         eigFreqOpt = 1; normOpt = 1; objOpt = 1;
    case 1.2
         eigFreqOpt = 2; normOpt = 1; objOpt = 1;
    case 2.0
         eigFreqOpt = 0; normOpt = 1; objOpt = 2;
    case 2.1
         eigFreqOpt = 1; normOpt = 1; objOpt = 2;
    case 2.2
         eigFreqOpt = 2; normOpt = 1; objOpt = 2;
    case 2.3
         eigFreqOpt = 0; normOpt = 2; objOpt = 2;
    case 2.4
         eigFreqOpt = 1; normOpt = 2; objOpt = 2;
    case 2.5
         eigFreqOpt = 2; normOpt = 2; objOpt = 2;
    otherwise
        error('\nWrong option for objective function!');
end



%% Structure parameter
N = size( structModel.K0, 1 );
n_alpha = length( structModel.K_j ) ;

% Identify the unmeasured DOFs
expModes.n_meas = length(expModes.measDOFs);
if (size(expModes.measDOFs, 1) < size(expModes.measDOFs, 2))
    expModes.measDOFs = expModes.measDOFs';
end
unmeasDOFs = setdiff( (1 : N)', expModes.measDOFs);
expModes.n_unmeas = length(unmeasDOFs);
expModes.n_modes = length(expModes.lambdaExp); % Number of measured modes

% Reordering stiffness, mass, influence matrices by measured and unmeasured
% DOFs
orderIndex = [expModes.measDOFs; unmeasDOFs];
structModel.M0 = structModel.M0(orderIndex, orderIndex);
structModel.K0 = structModel.K0(orderIndex, orderIndex);
for i = 1 : n_alpha
    structModel.K_j{i} = structModel.K_j{i}(orderIndex, orderIndex);
end


expModes.n_meas = length( expModes.measDOFs ); % # of measured DOFs
expModes.n_modes = length( expModes.lambdaExp ); % # of available exp-modes

% Add a new field K into structModel, which represents the stiffness matrix
% with current alpha values.
n_alpha = length(structModel.K_j);
structModel.K = structModel.K0;
for i = 1 : n_alpha
    structModel.K = structModel.K + x(i) * structModel.K_j{i};
end

if ( updatingOpts.formID ~= 3 )
    % When updatingOpts.formID equals 3, the modal dynamic residual
    % formulation does not require solving eigs. 
    numSimModes = max( updatingOpts.simModesForExpMatch );
    eigsOpts.tol = eps;
    [psi, lambda] = eigs( structModel.K, structModel.M0, numSimModes, ...
        'sm', eigsOpts );
    [lambda, ind] = sort( diag(lambda), 'ascend' );
    psi = psi(:, ind);
    % psi_m: entries in simulated mode shape psi that correspond to the
    % instrumented/measured DOFs.
    psi_m = psi(1 : expModes.n_meas, :);
    if ( updatingOpts.modeMatch == 1 )
        psiTemp = psi_m;
        matchedModeIndex = zeros(expModes.n_modes, 1);
        for i = 1 : expModes.n_modes
            [~, matchedModeIndex(i)] = max( MAC(expModes.psiExp(:, i), psiTemp) );
            % Wipe out the matched mode so that it is not again compared
            % with next experimental mode.
            psiTemp(:, matchedModeIndex(i)) = zeros(expModes.n_meas, 1);
        end
    elseif ( updatingOpts.modeMatch == 2 )
        matchedModeIndex = updatingOpts.simModesForExpMatch;
    else
        error ('\nWrong option for variable updatingOpts.modeMatch (1 or 2).');
    end
    
    % Create simModes variable to contain simulated modal properties that
    % will be used to match experimental modes when evaluating objective
    % function of Jacobians.
    %   lambda (n_modes x 1) - simulated eigenvalue
    %   psi_m  (n_meas x n_modes) - simulated mode shape vector at measured DOFs
    %   psi    (N x n_modes) - simulated mode shape vector at all DOFs
    simModes.psi_m = psi_m(:, matchedModeIndex);
    simModes.psi = psi(:, matchedModeIndex);
    simModes.lambda = lambda(matchedModeIndex);

    for i = 1 : expModes.n_modes
        if(~isempty(find(expModes.psiExp(:,i) == 1, 1)))
            expModes.q(i) = find(expModes.psiExp(:,i) == 1);
        else
            [~, expModes.q(i)] = max( abs( expModes.psiExp(:,i) ) );
        end
    end  
    
if( updatingOpts.formID < 2.3)
    % Normalize mode shape vector so that the maximum entry magnitude = 1.
    % The entry index is denoted as q(i).
    for i = 1 : expModes.n_modes
        if(~isempty(find(expModes.psiExp(:,i) == 1, 1)))
            expModes.q(i) = find(expModes.psiExp(:,i) == 1);
        else
            expModes.psiExp(:,i) = expModes.psiExp(:,i) / expModes.psiExp(expModes.q(i), i);
        end
    end
else
    % Normalize mode shape vector so that the length = 1.
    for i = 1 : expModes.n_modes
        expModes.psiExp(:,i) = expModes.psiExp(:,i) / norm(expModes.psiExp(:,i));
    end
end    

    if( updatingOpts.formID < 2.3)
        % Normalize mode shape vector so that the maximum entry magnitude = 1.
        % The entry index is denoted as q(i).
        for i = 1 : expModes.n_modes
            simModes.psi_m(:,i) = simModes.psi_m(:,i) / simModes.psi_m(expModes.q(i), i);
            simModes.psi(:,i) = simModes.psi(:,i) / simModes.psi(expModes.q(i), i);
        end
    else
        % Normalize mode shape vector so that the length = 1.
        for i = 1 : expModes.n_modes
            simModes.psi(:,i) = simModes.psi(:,i) / norm(simModes.psi_m(:, i));
            simModes.psi_m(:,i) = simModes.psi_m(:,i) / norm(simModes.psi_m(:, i));
        end
    end
    
end



%% Calculate Jacobian and Hessian matrix of residual
omegaSim = sqrt( simModes.lambda );
omegaExp = sqrt( expModes.lambdaExp );

n_meas = expModes.n_meas;
n_modes = expModes.n_modes;

modalMass = zeros( n_modes, 1 );
for i = 1 : n_modes
    modalMass(i) = simModes.psi(:, i)' * structModel.M0 * simModes.psi(:, i);
end

% dLambda will store values for d_lambda_i / d_alpha_j
dLambda = zeros(n_modes, n_alpha);
% ddLambda will store values for dd_lambda_i / d_alpha_j_d_alpha_k
ddLambda = zeros(n_alpha, n_alpha, n_modes);

% d_ri_eigFreqTerm will store the first part of d_r_i / d_alpha, i.e.
% the part involving the derivative of eigenvalue | angular
% frequency " ordinary frequency over alpha. For example, when eiegenvalue
% is used, the formulation is:
%       -weight_Lambda_i * D_alpha(Lambda_i) / Lambda_i^EXP
d_ri_eigFreqTerm = zeros( n_modes, n_alpha);
% dd_ri_eigFreqTerm will store the first part of dd_r_i / dd_alpha, i.e.
dd_ri_eigFreqTerm = zeros( n_alpha, n_alpha, n_modes);

% The 3rd dimension corresponds to i in Psi_i^m -- the mode index
dPsi_m_j = zeros(n_meas, n_alpha, n_modes);
dPsi_m_k = zeros(n_meas, n_alpha, n_modes);

if (objOpt == 1)
    jac_r = zeros(2 * n_modes, n_alpha);
    hes_r = zeros(n_alpha, n_alpha, n_modes * 2);
else
    if(normOpt == 1)
        jac_r = zeros(n_meas * n_modes, n_alpha);
        dd_ri_psi = zeros(n_alpha, n_alpha, (n_meas - 1) * n_modes);
        hes_r = zeros(n_alpha, n_alpha, n_modes * n_meas);
    else
        jac_r = zeros((n_meas + 1) * n_modes, n_alpha);
        dd_ri_psi = zeros(n_alpha, n_alpha, n_meas * n_modes);
        hes_r = zeros(n_alpha, n_alpha, n_modes * (n_meas + 1));
    end
end



for i = 1 : n_modes
    for j = 1 : n_alpha
        % eigenvalue first derivative
        dLambda(i,j) = simModes.psi(:,i)' * structModel.K_j{j} *...
            simModes.psi(:,i) / modalMass(i);
    end
end



for i = 1 : n_modes
    ddPsi_dAlpha_j_dAlpha_k = zeros(N,1);
    ddPsi_mr = zeros(n_alpha, n_alpha, n_meas - 1);
    ddPsi_m = zeros(n_alpha, n_alpha, n_meas);

    B = structModel.K - simModes.lambda(i) * structModel.M0;
    % The maximum entry of Psi_m is normalized to 1
    P_i = setdiff(1 : N, expModes.q(i));
    Q_i = P_i(1 : n_meas - 1);
    B = sparse(B(P_i, P_i));
    dcompB = decomposition(B); % factorization
    
    for j = 1 : n_alpha
        % eigenvector first derivative
        dPsi_dAlpha_j = zeros(N, 1);
        b_j = dLambda(i,j) * structModel.M0 * simModes.psi(:, i) -...
            structModel.K_j{j} * simModes.psi(:, i);
        % Cross out q_i-th row in B and b, and q_i-th column in B
        b_j = sparse(b_j(P_i));
        dPsi_dAlpha_j(P_i) = dcompB \ b_j;
        % The q_i-th entry of dPsi_m remains as 0 from initialization
        dPsi_m_j(Q_i, j, i) = dPsi_dAlpha_j(Q_i);
        
        for k = 1 : n_alpha
            dPsi_dAlpha_k = zeros(N, 1);
            b_k = dLambda(i,k) * structModel.M0 * simModes.psi(:, i) -...
                structModel.K_j{k} * simModes.psi(:, i);
            % Cross out q_i-th row in B and b, and q_i-th column in B
            b_k = b_k(P_i);
            dPsi_dAlpha_k(P_i) = dcompB \ sparse(b_k);
            % The q_i-th entry of dPsi_m remains as 0 from initialization           
            dPsi_m_k(Q_i, k, i) = dPsi_dAlpha_k(Q_i);
            
            % eigenvalue second derivative
            ddLambda(k,j,i) = ...
                (simModes.psi(:,i)' * (structModel.K_j{j} - dLambda(i,j) * structModel.M0 ) * dPsi_dAlpha_k ...
                + simModes.psi(:,i)' * (structModel.K_j{k} - dLambda(i,k) * structModel.M0 ) * dPsi_dAlpha_j) / modalMass(i);
            
            if eigFreqOpt == 0
                d_ri_eigFreqTerm(i,j) = - dLambda(i,j) * expModes.lambdaWeights(i) ...
                    / expModes.lambdaExp(i) ;
                dd_ri_eigFreqTerm(k,j,i) = - ddLambda(k,j,i) * expModes.lambdaWeights(i) ...
                    / expModes.lambdaExp(i);
                
            elseif eigFreqOpt == 1 || eigFreqOpt == 2
                d_ri_eigFreqTerm(i,j) = - dLambda(i,j) * expModes.lambdaWeights(i) ...
                    / (omegaExp(i) * 2 * omegaSim(i));
                dd_ri_eigFreqTerm(k,j,i) = - ddLambda(k,j,i) * expModes.lambdaWeights(i) ...
                    / (omegaExp(i) * 2 * omegaSim(i));
            end
            
            
            
            % eigenvector second derivative
            b_jk = ddLambda(k,j,i) * structModel.M0 * simModes.psi(:,i)...
                - (structModel.K_j{j} - dLambda(i,j) * structModel.M0) * dPsi_dAlpha_k ...
                - (structModel.K_j{k} - dLambda(i,k) * structModel.M0) * dPsi_dAlpha_j;
            b_jk = b_jk(P_i);
            ddPsi_dAlpha_j_dAlpha_k(P_i) = dcompB \ b_jk;
            
            for l = 1 : length(Q_i)
                ddPsi_mr(k, j, l) = ddPsi_dAlpha_j_dAlpha_k(Q_i(l));
            end
            ddPsi_m(k, j, Q_i) = ddPsi_dAlpha_j_dAlpha_k(Q_i);
            
        end
    end
    
    psiExp = expModes.psiExp(:,i);
    psiSim = simModes.psi_m(:,i);
    
    if (objOpt == 1)
        % For the MAC value formulation, d_ri_MAC is the second part of
        % d_r_i/d_alpha that involves MAC value:
        MACValue = MAC(psiExp, psiSim);
        d_ri_MAC = (-expModes.psiWeights(i) / sqrt(MACValue)) * ...
            (psiExp' / (psiExp' * psiSim) - psiSim' / norm(psiSim)^2) ...
            * dPsi_m_j(:, :, i);
        jac_r((i - 1) * 2 + 1 : i * 2, :) = [d_ri_eigFreqTerm(i,:); d_ri_MAC];
        
        dd_ri_MAC = zeros(n_alpha, n_alpha);
        for l = 1 : n_meas
            temp = (-expModes.psiWeights(i) / sqrt(MACValue)) * (psiExp' / (psiExp' * psiSim) - psiSim' / norm(psiSim)^2);
            dd_ri_MAC = dd_ri_MAC + temp(l) * ddPsi_m(:, :, l);
        end
        hes_r(:, :, (2 * i - 1)) = dd_ri_eigFreqTerm(:, :, i);
        hes_r(:, :, 2 * i) = dd_ri_MAC;
        
    elseif  (objOpt == 2)
        % For the eigenvector difference formulation, d_ri_psi is the
        % second part of d_r_i/d_alpha that involves eigenvectors:
        %      -Qi * D_alpha(Psi_i^m) * weight_Psi_i
        if(normOpt == 1)
            d_ri_psi = -dPsi_m_j(Q_i, :, i) * expModes.psiWeights(i);
            jac_r((i - 1) * n_meas + 1 : i * n_meas, :) = ...
                [d_ri_eigFreqTerm(i,:); d_ri_psi];
            
            dd_ri_psi = -ddPsi_mr * expModes.psiWeights(i);
            hes_r(:,:,i + (i - 1) * (n_meas - 1)) = [dd_ri_eigFreqTerm(:,:,i)];
            hes_r(:,:,1 + i + (i - 1) * (n_meas - 1) : i + i * (n_meas - 1)) = dd_ri_psi(:,:,:);
            
        elseif(normOpt == 2)
            d_ri_psi = -expModes.psiWeights(i) * ...
                ((psiExp' * psiSim* eye(n_meas) + psiSim *...
                psiExp') / norm(psiSim)^2 - 2 * (psiSim * psiExp' * (psiSim * psiSim'))...
                / norm(psiSim)^4) * dPsi_m_j(:, :, i);
            jac_r((i - 1) * (n_meas + 1) + 1 : i * (n_meas + 1),:) = ...
                [d_ri_eigFreqTerm(i,:); d_ri_psi];
            
            for l = 1 : n_meas
                temp = -expModes.psiWeights(i) * ...
                    ((psiExp' * psiSim* eye(n_meas) + psiSim * psiExp') / norm(psiSim)^2 ...
                    - 2 * (psiSim * psiExp' * (psiSim * psiSim')) / norm(psiSim)^4);
                dd_ri_psi(:,:,l + (i - 1) * n_meas) = temp(l) * ddPsi_m(:, :, l);
                hes_r(:,:,i + (i - 1) * n_meas) = [dd_ri_eigFreqTerm(:,:,i)];
                hes_r(:,:,l + i + (i - 1) * n_meas) = dd_ri_psi(:,:,l + (i - 1) * n_meas);
            end            
        end
    end
end



%% Assemble Jacobian and Hessian matrix for objective function
r =  ModelUpdatingObjective(x, structModel, expModes, simModes, updatingOpts);
temp = 0;
for i = 1:length(r)
    temp = temp + r(i) * hes_r(:,:,i);
end
hess = 2 * (jac_r' * jac_r + temp);
end
