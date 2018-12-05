function [r, jac] = LsqnonlinObjJac(x, structModel, expModes, updatingOpts,optToolBox)

% function [r, jac] = ModelUpdatingObjJac(x, structModel, expModes, updatingOpts,optToolBox)
%
% function [r] = ModelUpdatingObjJac(x, structModel, expModes, updatingOpts,optToolBox)
%
%   (c) Yang Wang, Xinjun Dong, Dan Li (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2018
%
% Revision: 1.0
%
% For implementation with different optimization algorithm, this function 
% calculates the objective residual vector r(x) and returns as the first 
% argument. When the function is called with two  output arguments:
%   [r, jac] = ModelUpdatingObjJac(alpha, structModel, expModes, updatingOpts)
% The function evaluates the user-provided analytical Jacobian matrix
% (d_r/d_x), and returns in the second argument.
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
%       modeMatch - Option for the matching method between simulated and
%       experimental modes (default: 1)
%           1: Match by the MAC value between the pair of simulated and
%             experimental mode shape vectors.
%           2: Strictly match the first designated simulated mode with the
%             first experimental mode (see simModesForExpMatch on how to
%             designate the simulated modes), the second designated
%             simulated mode with the second experimental mode, etc. Only
%             use this option when we are confident all the lowest few
%             experimental modes are captured, i.e. there is no missing
%             or unmeasured/undetected mode from the experimental data.
%       simModesForExpMatch - designate simulated modes obtained from FE
%       model for matching with experimental modes
%          If modeMatch = 1,
%               Set simModesForExpMatch as an integer representing the
%               number of simulated modes that will be compared with
%               experimental modes for similarity matching by MAC value.
%               The matched pair will be used for evaluating objective
%               function value. (default: min(n_modes x 2, N))
%          If modeMatch = 2,
%               Set simModesForExpMatch as a (n_modes x 1) array.  For
%               evaluating the objective function, the first experimental
%               mode will be matched with simModesForExpMatch(1)-th
%               simulated mode; the second experimental mode will be
%               matched with simModesForExpMatch(2)-th mode, etc.
%          x_ub - upper bounds of updating variables
%               formID < 3 - n_alpha x 1 (default: [])
%               formID = 3 - (n_alpha + n_unmeas x n_modes) x 1 (default: [])
%               Here n_unmeas refers to the number of unmeasured DOFs
%          x_lb - lowe bounds of updating variables
%               formID < 3 - n_alpha x 1 (default: [])
%               formID = 3 - (n_alpha + n_unmeas x n_modes) x 1 (default: [])
%
%          WARNING: when using Levenberg-Marquardt optimization algorithm in MATLAB,
%           setting the upper and lower bounds of updating variables has no
%           effect because the MATLAB L-M implementation does not accept
%           bounds. The optimization may provide an infeasible
%           out-of-the-bound solution. The user needs to verify the
%           feasibility of the solution.
%    optToolBox - option of optimization toolbox
%            'lsqnonlin'
%            'fmincon'
%            'Gauss-Newton'
%            'L_M'
%
% Output:
%	r: the objective residual vector r(x)
%   jac: the Jacobian matrix of the residual vector d_r / d_x

expModes.n_meas = length( expModes.measDOFs ); % # of measured DOFs
expModes.n_modes = length( expModes.lambdaExp ); % # of available exp-modes

% Create simModes variable to contain simulated modal properties that
% will be used to match experimental modes when evaluating objective
% function of Jacobians.
%   Lambda (n_modes x 1) - simulated eigenvalue
%   psi_m  (n_meas x n_modes) - simulated mode shape vector at measured DOFs
%   psi    (N x n_modes) - simulated mode shape vector at all DOFs
simModes = struct('psi_m',[],'psi',[],'Lambda',[]);
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
    
    simModes.psi_m = psi_m(:, matchedModeIndex);
    simModes.psi = psi(:, matchedModeIndex);
    simModes.Lambda = lambda(matchedModeIndex);
    
    if( updatingOpts.formID < 2.3 || updatingOpts.formID >= 4)
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

r =  ModelUpdatingObjective(x, structModel, expModes, simModes, updatingOpts);
r = sparse(r);

if nargout > 1
    jac = ModelUpdatingJacobian(x, structModel, expModes, simModes, updatingOpts);
    jac = sparse(jac);
end

% fmincon use scalar as objective function output
if(strcmp(optToolBox,'fmincon'))
	if nargout > 1
	    jac = 2 * jac' * r;
	end
    r = norm(r)^2;
end

end

