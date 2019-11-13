clear
close all
clc
addpath(genpath('StructuralDynamics'));

%% Strucutural parameters
n_modes = 1;                               % number of measured mode shapes
N = 4;                                     % DOF of the whole structure
masses = 12.060 / 386.088  * ones(N,1);    % mass value
iniSpring = 10 * ones(N,1);                % initial stiffness value   

%% Simulated actual structure
dmgLoc = 4;                                % damage location                              
alpha_act = [-0.1];                        % relative change from nominal 
                                           % stiffness value 
actSpring = iniSpring;
for i = 1:length(dmgLoc)
    actSpring(dmgLoc(i)) = actSpring(dmgLoc(i)) * (1 + alpha_act(i));
end

%% Actual strucuture
M0 = assemblyM(masses, N);
K_act = assemblyK(actSpring, N);
[psiExp,lambdaExp] = eigs(K_act, M0, 1,'sm') ;
[lambdaExp,dummyInd] = sort((diag(lambdaExp)),'ascend') ;
measDOFs = [1 2 3];
num_measDOFs = length(measDOFs);
unmeasDOFs = setdiff(1:N,measDOFs);
num_unmeasDOFs = length(unmeasDOFs);
lambdaExp = lambdaExp(1:n_modes);
modeIndex = 1 : n_modes;
psiExp = psiExp(:,dummyInd);

psiExp_m = psiExp(measDOFs,1:n_modes) ;     % simulated measurement
psiExp_u = psiExp(unmeasDOFs,1:n_modes) ;    
for i = 1:n_modes
    [~,q(i)] = max(abs(psiExp_m(:,i)));
    psiExp_u(:,i) = psiExp_u(:,i) / psiExp_m(q(i),i);
    psiExp_m(:,i) = psiExp_m(:,i) / psiExp_m(q(i),i);
end

%% Initial strucuture modes
M0 = assemblyM(masses, N);
K0 = assemblyK(iniSpring, N);

%% Assemble sensitivty matrix
i = 1 ;
K_j_Total(:,:,i) = zeros(N) ;
K_j_Total(i,i,i) = iniSpring(i);

for i = 2 : 4
    K_j_Total(:,:,i) = zeros(N) ;
    K_j_Total(i-1,i-1,i) = iniSpring(i);
    K_j_Total(i-1,i,i) = -iniSpring(i);
    K_j_Total(i,i-1,i) = -iniSpring(i);
    K_j_Total(i,i,i) = iniSpring(i);
end

K_j_ = K_j_Total(:,:,dmgLoc);


n_alpha = size(K_j_,3);


for i = 1 : n_alpha
    K_j{i} = sparse(K_j_(:,:,i));
end

%% Assemble structure matrices
structModel.M0 = M0;
structModel.K0 = K0;
structModel.K_j = K_j;

%% Optimization structure parameter;
optimzOpts.toolBox = 'fmincon';
optimzOpts.optAlgorithm = 'interior-point';
optimzOpts.gradSel = 'on';
optimzOpts.maxIter = 1e3;
optimzOpts.maxFunEvals = 3e5;


%% Model updating parameter
updatingOpts.formID = 3.0;       % 1: Modal property diff (MAC) ;
                                 % 2: Modal property diff (V_mDiff);
                                 % 3: Modal dynamic residual;
updatingOpts.modeMatch = 2;      % 1: Without forced matching;
                                 % 2: With forced matching;
updatingOpts.simModesForExpMatch = modeIndex;
if(updatingOpts.formID > 2.0 && updatingOpts.formID < 4.0)
    % Optimizaiton variable for modal dynamic residual formulation
    updatingOpts.x_lb = [-ones(n_alpha,1);-3 * ones(num_unmeasDOFs * n_modes,1)];
    updatingOpts.x_ub = [ones(n_alpha,1);  3 * ones(num_unmeasDOFs * n_modes,1)];
    
else
    % Optimizaiton variable for modal property difference formulations
    updatingOpts.x_lb = -ones(n_alpha,1);
    updatingOpts.x_ub =  ones(n_alpha,1);
end

%% Simulate "experimental data"
expModes.lambdaExp = lambdaExp;
expModes.psiExp = psiExp_m;
expModes.measDOFs = measDOFs;
expModes.lambdaWeights = ones(n_modes,1);
if(updatingOpts.formID == 2)
    expModes.psiWeights = ones(num_measDOFs,n_modes);
else
    expModes.psiWeights = ones(n_modes,1);
end

expModes.resWeights = ones(n_modes,1);

%% MultiStart optimization
numRuns = 10;
randSeed = 2;
filename = ['Sim_form' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_' optimzOpts.optAlgorithm '.mat'];

MultiRunModelUpdating




