function updtRslts = StructModelUpdating (structModel, expModes, ...
    updatingOpts, optimzOpts)
 
% function updtResults = StructModelUpdating (structModel, expModes, ...
%   updatingOpts, optimzOpts)
%
% function updtResults = StructModelUpdating (structModel, expModes, ...
%   [], optimzOpts)
%
% function updtResults = StructModelUpdating (structModel, expModes, ...
%   updatingOpts)
%
% function updtResults = StructModelUpdating (structModel, expModes)
%
%   Yang Wang, Xinjun Dong, Dan Li
%   School of Civil and Environmental Engineering
%   Georgia Institute of Technology
%   2018
%
% Revision: 1.0
%
% This function performs one run of the finite element model updating using
% frequency domain modal properties. The starting point for this run can be
% provided as optimzOpts.x0.
%
% Input:
%   structModel - a MATLAB structure array with following fields of
%   structural model information:
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
%
%   optimzOpts - optimization options. The current revision supports MATLAB
%          lsqnonlin and fmincon function.
%       maxIter - maximum iterations of optimization process (default: 400)
%       maxFunEvals - maximum number of function evaluations
%           allowed (default: 100 x n_alpha)
%       tolFun - termination tolerance on the value change of the objective
%           function between two iterations (default: 1e-6)
%       tolX - termination tolerance on the value change of the optimization
%           vector variables between two iterations (default: 1e-6)
%       gradSel - selection for gradient calculation
%          'on': calculate search gradient with user-defined Jacobian
%             matrix. With r representing the residual vector whose length
%             square is the objective function value, The gradient
%                 [d_f / d_alpha]' = [(d_f / d_r) * ((d_r / d_alpha)]'
%          'off': let MATLAB numerically calculate gradient matrix
%               using finite difference method (default)
%       toolBox - Optimization toolbox
%            'lsqnonlin' (default)
%            'fmincon'
%       optAlgorithm - optimization algorithm
%           lsqnonlin
%            'trust-region-reflective' algorithm
%            'Levenberg-Marquardt' algorithm (default)
%           fmincon
%            'trust-region-reflective' algorithm
%            'interior-point' algorithm
%       x0 - initial value for updating variables
%            formID < 3 - n_alpha x 1
%            formID = 3 - (n_alpha + n_unmeas x n_modes) x 1
%            (default: zero vector)
%%
% Output:
%   A structure array with following fields:
%      xOpt - optimal value of updating variables
%      fvalOpt - optimal value of objective function value
%      exitFlag - exit flag of MATLAB lsqnonlin (check MATLAB lsqnolin help
%         for detail: https://www.mathworks.com/help/optim/ug/lsqnonlin.html)
%      gradient - gradient of objective function at alphaOpt

% Structure Parameter
N = size(structModel.M0);
n_alpha = length(structModel.K_j);

% Identify the unmeasured DOFs
expModes.n_meas = length(expModes.measDOFs);
if (size(expModes.measDOFs, 1) < size(expModes.measDOFs, 2))
    expModes.measDOFs = expModes.measDOFs';
end
unmeasDOFs = setdiff( (1 : N)', expModes.measDOFs);
expModes.n_unmeas = length(unmeasDOFs);
expModes.n_modes = length(expModes.lambdaExp); % Number of measured modes

if(nargin < 3 || isempty(updatingOpts))
    % Default model updating options if not provided
    updatingOpts = struct('formID', 1, 'modeMatch', 1,...
        'simModesForExpMatch', min([expModes.n_modes * 2, N]),...
        'x_lb', [], 'x_ub',[]);
else
    % Default model updating options for missing fields
    if (~isfield(updatingOpts, 'formID'))
        updatingOpts.formID = 1;
    end
    
    if (~isfield(updatingOpts, 'modeMatch'))
        updatingOpts.modeMatch = 1;
    end
    
    if (~isfield(updatingOpts, 'simModesForExpMatch'))
        updatingOpts.simModesForExpMatch = min([expModes.n_modes * 2, N]);
    end
    
    if (updatingOpts.modeMatch == 1)
        if (length(updatingOpts.simModesForExpMatch) ~= 1) || ...
                (updatingOpts.simModesForExpMatch < expModes.n_modes)
            error(['Please provide correct number for updatingOpts.simModesForExpMatch. ' ...
                'When modeMatch == 1, simModesForExpMatch should be a scalar representing '...
                'the number of simulated modes for experimental matching. The scalar value ' ...
                'should be no smaller than n_modes = %d.'], expModes.n_modes);
        end
    elseif (updatingOpts.modeMatch == 2 && length(updatingOpts.simModesForExpMatch) < expModes.n_modes)
        error(['Please provide correct number of updatingOpts.simModesForExpMatch if modeMatch is set to be 2.' ...
            '\nThe number of simulated modes for experimental matching should be the same as expModes.n_modes = %d.'],...
            expModes.n_modes);
    end
    
    if (~isfield(updatingOpts, 'x_lb'))
        updatingOpts.x_lb = [];
    end
    
    if(~isfield(updatingOpts, 'x_ub'))
        updatingOpts.x_ub = [];
    end
end

if(updatingOpts.formID < 3)
    n_x = n_alpha;
    
else
    n_x = n_alpha + expModes.n_unmeas * expModes.n_modes;
end


if(nargin < 4)
    % Default optimization options if not provided
    optimzOpts = struct ('maxIter', 400, 'maxFunEvals', 100 * n_x, ...
        'tolFun', 1e-6, 'tolX', 1e-6, 'gradSel', 'off',...
        'optAlgorithm', 'Levenberg-Marquardt', 'x0', zeros(n_x, 1));
else
    % Default optimization options
    if (~isfield(optimzOpts, 'maxIter'))
        optimzOpts.maxIter = 400;
    end
    
    if (~isfield(optimzOpts, 'maxFunEvals'))
        optimzOpts.maxFunEvals = 100 * n_x;
    end
    
    if (~isfield(optimzOpts, 'tolFun'))
        optimzOpts.tolFun = 1e-6;
    end
    
    if (~isfield(optimzOpts, 'tolX'))
        optimzOpts.tolX = 1e-6;
    end
    
    if (~isfield(optimzOpts, 'gradSel'))
        optimzOpts.gradSel = 'off';
    end
    
    if (~isfield(optimzOpts, 'optAlgorithm'))
        optimzOpts.optAlgorithm = 'Levenberg-Marquardt';
    end
    
     if (~isfield(optimzOpts, 'toolBox'))
        optimzOpts.toolBox = 'lsqnonlin';
    end
    
        
    if (~isfield(optimzOpts, 'x0'))
        optimzOpts.x0 = zeros(n_x, 1);
    end
    
    if (length(optimzOpts.x0) ~= n_alpha && updatingOpts.formID < 3)
        error(['\nThe number of sensitivity matrices structModel.K_j does not ' ...
            'match with the number of updating variables n_alpha = %d'], n_alpha);
    end
end

% Reordering stiffness, mass, influence matrices by measured and unmeasured
% DOFs
orderIndex = [expModes.measDOFs; unmeasDOFs];
structModel.M0 = structModel.M0(orderIndex, orderIndex);
structModel.K0 = structModel.K0(orderIndex, orderIndex);
for i = 1 : n_alpha
    structModel.K_j{i} = structModel.K_j{i}(orderIndex, orderIndex);
end

if( updatingOpts.formID < 2.3)
    % Normalize mode shape vector so that the maximum entry magnitude = 1.
    % The entry index is denoted as q(i).
    for i = 1 : expModes.n_modes
        if(~isempty(find(expModes.psiExp(:,i) == 1, 1)))
            expModes.q(i) = find(expModes.psiExp(:,i) == 1);
        else
            [~, expModes.q(i)] = max( abs( expModes.psiExp(:,i) ) );
            expModes.psiExp(:,i) = expModes.psiExp(:,i) / expModes.psiExp(expModes.q(i), i);
        end
    end
else
    % Normalize mode shape vector so that the length = 1.
    for i = 1 : expModes.n_modes
        expModes.psiExp(:,i) = expModes.psiExp(:,i) / norm(expModes.psiExp(:,i));
    end
end

if( (updatingOpts.formID >= 2 && updatingOpts.formID < 3) || updatingOpts.formID >= 5) 
    if(size(expModes.psiWeights,1) <  expModes.n_meas)
        psiWeights = expModes.psiWeights';
        expModes.psiWeights = zeros(expModes.n_meas, expModes.n_modes);
        for i = 1:expModes.n_meas
            expModes.psiWeights(i,:) =  psiWeights .* ones(1, expModes.n_modes);
        end
    end
    
end
% Accommodate different variable names in different MATLAB versions.
vers_temp = version( '-release' );
MAT_version = str2num( vers_temp(1 : 4) );
if (MAT_version >= 2016)
    maxIterName = 'MaxIterations';
    maxFunName = 'MaxFunctionEvaluations';
else
    maxIterName = 'maxIter';
    maxFunName = 'maxFunEvals';
end

% If the bounds are not empty, MATLAB forces to use
% trust-region-reflective algorithm.
if strcmp(optimzOpts.optAlgorithm, 'Levenberg-Marquardt')
    updatingOpts.x_lb = [];
    updatingOpts.x_ub = [];
end

fun = @(x) OptmzObjJac(x, structModel, expModes, updatingOpts,optimzOpts.toolBox);

if(strcmp(optimzOpts.toolBox,'lsqnonlin'))
    options = optimoptions( 'lsqnonlin', 'tolFun', optimzOpts.tolFun, 'tolX', optimzOpts.tolX,...
        'Algorithm', optimzOpts.optAlgorithm, 'Disp', 'iter', 'Jacobian', optimzOpts.gradSel,...
        maxIterName, optimzOpts.maxIter, maxFunName, optimzOpts.maxFunEvals );
    if strcmp(optimzOpts.optAlgorithm, 'Levenberg-Marquardt')
        updatingOpts.x_lb = [];
        updatingOpts.x_ub = [];
    end
    % If the bounds are not empty, MATLAB forces to use
    % trust-region-reflective algorithm.
    [updtRslts.xOpt, updtRslts.fvalOpt, updtRslts.residual, updtRslts.exitFlag, ...
        updtRslts.output, ~, Jac_temp] = lsqnonlin( fun, optimzOpts.x0, ...
        updatingOpts.x_lb, updatingOpts.x_ub, options );
    updtRslts.gradient = 2 * full(Jac_temp)' * updtRslts.residual;
elseif(strcmp(optimzOpts.toolBox,'fmincon'))
    if(strcmp(optimzOpts.gradSel,'on'))
        optimzOpts.gradSel = true;
    else
        optimzOpts.gradSel = false;
    end
    options = optimoptions( 'fmincon', 'tolFun', optimzOpts.tolFun, 'tolX', optimzOpts.tolX,...
        'Algorithm', optimzOpts.optAlgorithm, 'Disp', 'iter', 'SpecifyObjectiveGradient',optimzOpts.gradSel,...
        maxIterName, optimzOpts.maxIter, maxFunName, optimzOpts.maxFunEvals );
    [updtRslts.xOpt, updtRslts.fvalOpt, updtRslts.exitFlag, updtRslts.output,~,updtRslts.gradient] = fmincon( fun, optimzOpts.x0,[],[],[],[], ...
        updatingOpts.x_lb, updatingOpts.x_ub,[], options );
end
