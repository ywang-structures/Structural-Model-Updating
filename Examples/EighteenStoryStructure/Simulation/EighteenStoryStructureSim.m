clear
close all
warning('off')

%% Basic parameters
N = 18; % # of DOFs of the whole structure
n_modes = 4; % # of measured experimental modes
measDOFs = [1; 4; 7; 9; 12; 15; 18]; % Measured DOFs
modeIndex = 1:n_modes; % Indexes of simulated modes for matching exp. modes 

%% Assemble structure matrices
stiffnesses = [1155; 1092; 1073; 1028; 1028; 990;
                963;  938;  876;  840;  824; 788;
                712;  660;  619;  562;  491; 363]*10^5; % N/m
weights = [208; 208; 208; 208; 208; 208;
           208; 208; 208; 208; 208; 206;
           206; 206; 206; 206; 206; 202]*10^3; % N
masses = weights/9.8; % kg        
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

%% Simulate "experimental data"
% actual values for stiffness updating variables, each alpha represents
% relative change from nominal stiffness parameter.
alpha_act = [0.05;  0.05; -0.05; -0.10;  0.10; -0.15;
             0.15;  0.25; -0.10;  0.20;  0.30;  0.25;
            -0.15;  0.05; -0.15;  0.10;  0.20;  0.20;];
K_act = sparse(zeros(N, N));
for i = 1:n_alpha
    K_act = K_act + K_j{i} * (1 + alpha_act(i)); % actual structure
end        
[psiExpAll,lambdaExp] = eigs(K_act,M0,n_modes,'sm');
[lambdaExp,dummyInd] = sort((diag(lambdaExp)),'ascend') ;
lambdaExp = lambdaExp(modeIndex);
psiExpAll = psiExpAll(:,dummyInd(modeIndex));
psi_m = psiExpAll(measDOFs,:);

% Normalize the mode shape vectors by maximum entry
for i = 1:n_modes
    [~,index] = max(abs(psi_m(:,i)));
    psi_m(:,i) = psi_m(:,i) / psi_m(index,i);
end

expModes.lambdaExp = lambdaExp;
expModes.psiExp = psi_m;
expModes.measDOFs = measDOFs;
unmeasDOFs = setdiff(1 : N, measDOFs);
num_measDOFs = length(measDOFs);
num_unmeasDOFs = length(unmeasDOFs);

%% Optimization settings;
optimzOpts.toolBox  = 'lsqnonlin'; % Optimization solver
optimzOpts.optAlgorithm = 'Levenberg-Marquardt'; % Optimization algorithm
optimzOpts.gradSel = 'on'; % on = anlaytical Jacobian; off = numerical Jacobian

%% Model updating settings
updatingOpts.formID = 1;         % 1: Modal property diff. formulation with MAC values
                                 % 2: Modal property diff. formulation with eigenvectors
                                 % 3: Modal dynamic residual formulation                               
updatingOpts.modeMatch = 2;      % 1: Without forced matching;
                                 % 2: With forced matching;
updatingOpts.simModesForExpMatch = modeIndex;

if(updatingOpts.formID < 3)
    % Optimization variable bounds for modal property difference formulations
    updatingOpts.x_lb = -ones(n_alpha,1);
    updatingOpts.x_ub =  ones(n_alpha,1);
else
    % Optimization variable bounds for modal dynamic residual formulation
    updatingOpts.x_lb = [-1.0*ones(n_alpha,1); -2 * ones(num_unmeasDOFs * n_modes,1)];
    updatingOpts.x_ub =  [1.0*ones(n_alpha,1);2 * ones(num_unmeasDOFs * n_modes,1)];    
end

% Weighting factors
expModes.lambdaWeights = ones(n_modes,1);
expModes.psiWeights = ones(n_modes,1);
expModes.resWeights = ones(n_modes,1);

%% MultiStart optimization
numRuns = 100;
randSeed = 1;
filename = ['EighteenStoryStructureSim_form' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_' optimzOpts.optAlgorithm '.mat'];

MultiRunModelUpdating
