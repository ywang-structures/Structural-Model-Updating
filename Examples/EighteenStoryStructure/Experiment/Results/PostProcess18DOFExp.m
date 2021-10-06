clc;clear;close all
font = 10;
figSize1 = [100, 100, 400, 220];
figSize2 = [100, 100, 230, 320];
set(0,'defaultAxesFontSize',font)
load('..\ExpModalInfo_18DOF.mat')

%% Load file
filename = 'EighteenStoryStructureExp_form2_JACon_Interior-point'
load(filename)
n_modes = 2; % number of modes used in model updating
[fvalmin,idxMin] = min(fval); % obtain the minimum obj value
alpha = x(1:18,idxMin); % obtain the corresponding stiffness updating variable

%% Initial model
N = 18;
k = [1155; 1092; 1073; 1028; 1028; 990;
                963;  938;  876;  840;  824; 788;
                712;  660;  619;  562;  491; 363]*10^2; % kN/m
weights = [208; 208; 208; 208; 208; 208;
           208; 208; 208; 208; 208; 206;
           206; 206; 206; 206; 206; 202]; % kN
masses = weights/9.8; % ton        
M0 = diag(masses);
K0 = makeK(k,N);
[psiSim, lambdaSim] = eig(K0, M0);
[lambdaSim,idxDummy] = sort((diag(lambdaSim)),'ascend') ;
lambdaSim = lambdaSim(1:n_modes);
freqSim = sqrt(lambdaSim)/2/pi;
psiSim = psiSim(:,idxDummy);
for i = 1:n_modes
    [~,q(i)] = max(abs(psiExp_m(:,i)));
    temp = psiSim(measDOFs(q(i)),i);
    psiSim(:,i) = psiSim(:,i) /temp;
end

%% Updated model
kupd = k .* (1 + alpha);
Kupd = makeK(kupd,N);
[psiUpd, lambdaUpd] = eig(Kupd, M0);
[lambdaUpd,idxDummy] = sort((diag(lambdaUpd)),'ascend') ;
lambdaUpd = lambdaUpd(1:n_modes);
freqUpd = sqrt(lambdaUpd)/2/pi;
psiUpd = psiUpd(:,idxDummy);
for i = 1:n_modes
    [~,q(i)] = max(abs(psiExp_m(:,i)));
    temp = psiUpd(measDOFs(q(i)),i);
    psiUpd(:,i) = psiUpd(:,i) /temp;
end

%% Comparison between the experimental structure, the initial model, and the updated model
freqExp = round(freqExp,3);
freqSim = round(freqSim,3);
freqUpd = round(freqUpd,3);
freqDiffSim = (freqSim - freqExp(1:n_modes))./freqExp(1:n_modes)*100;
freqDiffUpd = (freqUpd - freqExp(1:n_modes))./freqExp(1:n_modes)*100;

freqDiffSim % frequency diff (%) between the initial model and the experiemnt
freqDiffUpd % frequency diff (%) between the updated model and the experiemnt

%% Plot
figure('position', figSize1)
plot(1:length(fval),round(fval,2),'kx','MarkerSize',5);hold on
xlabel('Search cases')
ylabel({'Objective function value'})

figure('position', figSize2)
plot([0; real(psiExp_m(:,1))],[0;measDOFs'],'k*');hold on
plot([0; real(psiSim(:,1))],[0:N],'r--')
plot([0; real(psiUpd(:,1))],[0:N],'b')
xlabel('Eigenvector 1st mode')
ylabel('DOF')
legend('Experiment','Initial','Updated')
xlim([-1.5 1.5])
grid on

figure('position', figSize2)
plot([0; real(psiExp_m(:,2))],[0;measDOFs'],'k*');hold on
plot([0; real(psiSim(:,2))],[0:N],'r--')
plot([0; real(psiUpd(:,2))],[0:N],'b')
xlabel('Eigenvector 2nd mode')
ylabel('DOF')
legend('Experiment','Initial','Updated')
xlim([-1.5 1.5])
grid on

figure('position', figSize2)
h = barh(1:N,[k kupd],0.5);
h(1).FaceColor = 'r';
h(2).FaceColor = 'b';
xlabel({'Inter-story stiffness', '(kN/m)'})
ylabel('DOF')
ylim([0 19])
yticks([2:2:18])
xlim([0 max(kupd)])
legend('Initial','Updated')

figure('position', figSize2)
barh(1:N,alpha,0.4,'k');hold on
xlabel('Stiffness updating variable \alpha')
ylabel('DOF')
ylim([0 19])
xlim([-0.35 0.35])
xticks([-0.3:0.15:0.3])
yticks([1:1:18])
grid on
