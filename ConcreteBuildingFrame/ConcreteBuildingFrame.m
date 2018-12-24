clc;clear
close all
warning('off')


%% Actual values for stiffness updating variables, each alpha represents 
% relative change from nominal stiffness parameter.
alpha_act = [-0.1;  0.2;   0.2; -0.05; 0.2;  0.15;
             0.15; 0.10; -0.10; -0.15; 0.20; 0.15;];

%% Load Sensitivity matrix
LoadStructure
n_modes = 3; % Number of measured modes
modeIndex = 1 : n_modes; % Indexes of these measured modes

%% Assemble structure matrices
structModel.M0 = M0;
structModel.K0 = K0;
structModel.K_j = K_j;

%% Optimization structure parameter;
optimzOpts.tolFun = 1e-10;
optimzOpts.tolX = 1e-10;
optimzOpts.toolBox  = 'lsqnonlin';
optimzOpts.optAlgorithm = 'Levenberg-Marquardt';
optimzOpts.gradSel = 'on';
optimzOpts.maxIter = 1e3;
optimzOpts.maxFunEvals = 12e4;


%% Simulate "experimental data"
[psiExpAll, lambdaExp] = eigs(K_act, M0, n_modes, 'sm');

[lambdaExp,dummyInd] = sort(diag(lambdaExp), 'ascend') ;
lambdaExp = lambdaExp(modeIndex);
psiExpAll = psiExpAll(:, dummyInd(modeIndex));
psi_m = psiExpAll(measDOFs,:);

% Normalize the mode shape vectors by maximum entry
for i = 1: n_modes
    [~,index] = max(abs(psi_m(:,i)));
    psi_m(:,i) = psi_m(:,i) / psi_m(index,i);
end

expModes.lambdaExp = lambdaExp;
expModes.psiExp = psi_m;
expModes.measDOFs = measDOFs;
expModes.lambdaWeights = ones(n_modes,1);
expModes.psiWeights = ones(n_modes,1);


%% Model updating parameter
updatingOpts.formID = 3.0;       % 1: Modal property diff (MAC) ;
                                 % 2: Modal property diff (V_mDiff);
updatingOpts.modeMatch = 2;      % 1: Without forced matching;
                                 % 2: With forced matching;
updatingOpts.simModesForExpMatch = modeIndex;
updatingOpts.x_lb = -ones(n_alpha,1);
updatingOpts.x_ub =  ones(n_alpha,1);

%% MultiStart optimization
numRuns = 100;
randSeed = 3;
filename = ['ConcBuildFrm_form' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_' optimzOpts.optAlgorithm '.mat'];

MultiRunModelUpdating


