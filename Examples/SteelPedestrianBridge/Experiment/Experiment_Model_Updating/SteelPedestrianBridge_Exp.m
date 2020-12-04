%% 
clear
close all

LoadStructure;

modeIndex = [1 2 3 4 5];
n_modes = length(modeIndex);
%% Assemble structure parameter:
structModel.M0 = M0;
structModel.K0 = K0;
structModel.K_j = K_j;

%% Optimization structure parameter;
optimzOpts.tolFun = 1e-10;
optimzOpts.tolX = 1e-10;
optimzOpts.tolGrad = 1e-10;
% optimzOpts.toolBox = 'lsqnonlin';
% optimzOpts.optAlgorithm = 'trust-region-reflective';
optimzOpts.toolBox = 'fmincon';
optimzOpts.optAlgorithm = 'interior-point';

optimzOpts.gradSel = 'on';
optimzOpts.maxIter = 1e3;
optimzOpts.maxFunEvals = 3e3;

%% Model updating parameter
updatingOpts.formID = 2;       % 1: Modal property diff (MAC) ;
                               % 2: Modal property diff (V_mDiff);
                               % 3: Modal dynamic residual;
                               
updatingOpts.modeMatch = 2;    % 1: Without forced matching;
                               % 2: With forced matching;
                               
updatingOpts.simModesForExpMatch = modeIndex;

%% Experimental modal properties
load ExpModeInfo_FEM

expModes.lambdaExp = lambdaExp;
expModes.psiExp = psiExp_m;
expModes.measDOFs = measDOFs;

num_measDOFs = length(measDOFs);
num_unmeasDOFs = N - num_measDOFs;
unmeasDOFs = setdiff(1 : N, measDOFs);


if(updatingOpts.formID < 3)
    
    updatingOpts.x_lb = [-0.1 * ones(2, 1); -0.2 * ones(11,1);   -ones(4,1)];
    updatingOpts.x_ub = [ 0.1 * ones(2, 1);  0.4 * ones(11,1);  10 * ones(4,1) ];
    
else
    
    updatingOpts.x_lb = [-0.1 * ones(2, 1); -0.2 * ones(11,1);      -ones(4,1); -5 * ones(num_unmeasDOFs * n_modes,1)];
    updatingOpts.x_ub = [ 0.1 * ones(2, 1);  0.4 * ones(11,1);  10 * ones(4,1);  5 * ones(num_unmeasDOFs * n_modes,1) ];
end


for i = 1:length(modeIndex)
    if(updatingOpts.formID == 1 || updatingOpts.formID == 3)
        q = find(psiExp_m(:,i) == 1);
        indx = setdiff(1 : num_measDOFs, q);
        expModes.psiWeights(i,1) = 1/ norm(std_psiExp(indx,i));
    else
        expModes.psiWeights(:,i) = 1 ./ std_psiExp(:,i) / num_measDOFs;
    end
end

expModes.lambdaWeights = lambdaExp ./ std_lambdaExp(modeIndex);

% Weighitng factor for modal dynamic residual formualation
if(updatingOpts.formID == 3)
    %     expModes.resWeights =  weightLambda .* weightPsi;
    expModes.resWeights =  expModes.lambdaWeights .*  expModes.psiWeights;
end


start_Ch = 'SteelPedBridg_form';

filename = [start_Ch num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_' optimzOpts.optAlgorithm '.mat'];

% randSeed = 2;
numRuns = 100;
MultiRunModelUpdating
% 


