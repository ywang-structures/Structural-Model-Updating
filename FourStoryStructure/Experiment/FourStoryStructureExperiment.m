clear
close all
clc

%% Basic parameter
N = 4;                                    % DOF of the whole structure
n_modes = 2;                              % number of measured mode shapes
modeIndex = 1:n_modes;                    % Indexes of measured modes
measDOFs = [1;2;3];                       % Indexes of measured DOFs
unmeasDOFs = setdiff(1 : N,measDOFs)';    % Indexes of unmeasured DOFs
num_measDOFs = length(measDOFs);
num_unmeasDOFs = length(unmeasDOFs);

%% Load experimental data
load FourStoryStructure
lambdaExp = (test_frequency * 2 * pi).^2;
lambdaExp = lambdaExp(1 : n_modes);
psiExp_m = test_mode(measDOFs,1:n_modes);

for i = 1:n_modes
    % normlize q_i the maximum entry equal to 1
    [~,q(i)] = max(abs(psiExp_m(:,i)));
    psiExp_m(:,i) = psiExp_m(:,i) /  psiExp_m(q(i),i);
    % residual measured DOF other than the qth one
    measDOFsR(:,i) = setdiff(1:length(measDOFs),q(i),'stable')'; 
    psiExp_mR(:,i) = psiExp_m(measDOFsR(:,i),i);
end

%% Influence Matrix
stiffnesses = 0.01  * 1000 * ones(N,1);
masses = 12.060 / 386.088  * ones(N,1); 
i = 1 ;
K_j_(:,:,i) = zeros(N) ;
K_j_(i,i,i) = stiffnesses(i);

for i = 2 : N
    K_j_(:,:,i) = zeros(N) ;
    K_j_(i-1,i-1,i) = stiffnesses(i);
    K_j_(i-1,i,i) = -stiffnesses(i);
    K_j_(i,i-1,i) = -stiffnesses(i);
    K_j_(i,i,i) = stiffnesses(i);
end

n_alpha = size(K_j_,3);

for i = 1 : n_alpha
    K_j{i} = sparse(K_j_(:,:,i));
end
alpha_lb = -ones(n_alpha,1);
alpha_ub = -alpha_lb;

%% Initial Structure
rordIdx = [measDOFs;unmeasDOFs];
K0 = assemblyK(stiffnesses, N);
K0 = K0(rordIdx,rordIdx);
M0 = assemblyM(masses, N);
M0 = M0(rordIdx,rordIdx);

for i = 1 : n_alpha
    K_j{i} = K_j{i}(rordIdx,rordIdx);
end

[psiSimAll,lambdaSim] = eigs(K0,M0,n_modes,'sm');
[lambdaSim,dummyInd] = sort((diag(lambdaSim)),'ascend') ;
lambdaSim = lambdaSim(modeIndex);
psiSimAll = psiSimAll(:,dummyInd(modeIndex));
psiSim_m = psiSimAll(measDOFs,:);
psiSim_u = psiSimAll(unmeasDOFs,:);

for i = 1:n_modes
    temp = psiSim_m(q(i),i);
    psiSim_m(:,i) = psiSim_m(:,i) /temp;
    psiSim_u(:,i) = psiSim_u(:,i) / temp;
    psiSim_mR(:,i) = psiSim_m(measDOFsR(:,i),i);
end

%% Assemble structure matrices
structModel.M0 = M0;
structModel.K0 = K0;
structModel.K_j = K_j;

%% Optimization structure parameter;
optimzOpts.tolFun = eps^2;
optimzOpts.tolX = eps^2;
optimzOpts.toolBox = 'fmincon';
optimzOpts.optAlgorithm = 'interior-point';
optimzOpts.gradSel = 'on';
optimzOpts.maxIter = 1e3;
optimzOpts.maxFunEvals = 3e5;


%% Model updating parameter
updatingOpts.formID = 3.0;       % 1: Modal property diff (MAC) ;
                                 % 2: Modal property diff (V_mDiff);
updatingOpts.modeMatch = 2;      % 1: Without forced matching;
                                 % 2: With forced matching;
updatingOpts.simModesForExpMatch = modeIndex;
if(updatingOpts.formID > 2.0 && updatingOpts.formID < 4.0)
    % Optimizaiton variable for modal dynamic residual formulation
    updatingOpts.x_lb = [-ones(n_alpha,1);-2 * ones(num_unmeasDOFs * n_modes,1)];
    updatingOpts.x_ub = [ones(n_alpha,1);  2 * ones(num_unmeasDOFs * n_modes,1)];
    
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
numRuns = 1000;
randSeed = 2;
filename = ['Exp_form' num2str(updatingOpts.formID) '_Jac' optimzOpts.gradSel '_' optimzOpts.optAlgorithm '.mat'];

MultiRunModelUpdating




