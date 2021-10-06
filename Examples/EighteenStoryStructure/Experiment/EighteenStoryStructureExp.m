clear
close all
warning('off')

%% Basic parameter
N = 18; % # of DOFs of the whole structure
n_modes = 2; % # of measured experimental modes
modeIndex = 1:n_modes; % Indexes of simulated modes for matching exp. modes 

%% Assemble structure matrices
stiffnesses = [1155; 1092; 1073; 1028; 1028; 990;
                963;  938;  876;  840;  824; 788;
                712;  660;  619;  562;  491; 363]*10^2; % kN/m
weights = [208; 208; 208; 208; 208; 208;
           208; 208; 208; 208; 208; 206;
           206; 206; 206; 206; 206; 202]; % kN
masses = weights/9.8; % ton        
M0 = diag(masses);

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

K0 = sparse(zeros(N, N));
for i = 1:n_alpha
    K0 = K0 + K_j{i} ; % initial structure
end

structModel.M0 = M0;
structModel.K0 = K0;
structModel.K_j = K_j;

%% Load system identification data
load ExpModalInfo_18DOF
lambdaExp = lambdaExp(modeIndex);
psiExp_m = psiExp_m(:, modeIndex);
for i = 1:n_modes
    % normlize q_i the maximum entry equal to 1
    [~,q(i)] = max(abs(psiExp_m(:,i)));
    psiExp_m(:,i) = psiExp_m(:,i) /  psiExp_m(q(i),i);
end

expModes.lambdaExp = lambdaExp;
expModes.psiExp = psiExp_m;
expModes.measDOFs = measDOFs;
unmeasDOFs = setdiff(1 : N, measDOFs);
num_measDOFs = length(measDOFs);
num_unmeasDOFs = length(unmeasDOFs);

%% Optimization settings;
optimzOpts.toolBox  = 'fmincon'; % Optimization solver
optimzOpts.optAlgorithm = 'Interior-point'; % Optimization algorithm
optimzOpts.gradSel = 'on'; % on = anlaytical Jacobian; off = numerical Jacobian

%% Model updating settings
updatingOpts.formID = 2;         % 1: Modal property diff. formulation with MAC values
                                 % 2: Modal property diff. formulation with eigenvectors
                                 % 3: Modal dynamic residual formulation                               
updatingOpts.modeMatch = 2;      % 1: Without forced matching;
                                 % 2: With forced matching;
updatingOpts.simModesForExpMatch = modeIndex;

if(updatingOpts.formID < 3)
    % Optimization variable bounds for modal property difference formulations
    updatingOpts.x_lb = -0.3*ones(n_alpha,1);
    updatingOpts.x_ub =  0.3*ones(n_alpha,1);
else
    % Optimization variable bounds for modal dynamic residual formulation
    updatingOpts.x_lb = [-0.3*ones(n_alpha,1); -2 * ones(num_unmeasDOFs * n_modes,1)];
    updatingOpts.x_ub =  [0.3*ones(n_alpha,1);2 * ones(num_unmeasDOFs * n_modes,1)];    
end

expModes.lambdaWeights = ones(n_modes,1)*100; % assign 100 on eigenvalue difference
expModes.psiWeights = ones(n_modes,1);
expModes.resWeights = ones(n_modes,1);

%% MultiStart optimization
numRuns = 1000;
randSeed = 2;
filename = ['EighteenStoryStructureExp_form' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_' optimzOpts.optAlgorithm '.mat'];

MultiRunModelUpdating
