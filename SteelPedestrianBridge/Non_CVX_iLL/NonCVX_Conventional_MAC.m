% Updating for NEES: frequency domain
clear
close all

alpha_act = [0.05; 0.05; -0.05; -0.10; 0.10; -0.15;
             0.15; -0.05; -0.10; 0.10; -0.20;
             -0.30;0.60;-0.30;0.60;];
cd ..
LoadStructure;

n_alpha = length(alpha_act);
n_modes = 3;
modeIndex = 1:n_modes;
cd Non_CVX_iLL

%% Assemble structure paramter
structModel.M0 = M0;
structModel.K0 = K0;
structModel.K_j = K_j;

%% Simulate "experimental data"

[psiExp,lambdaExp] = eigs(K_act,M0,20,'sm');
[lambdaExp,dummyInd] = sort((diag(lambdaExp)),'ascend') ;

lambdaExp = lambdaExp(modeIndex);

psiExp = psiExp(:,dummyInd(modeIndex));

psiExp_m = psiExp(measDOFs,:);


for i = 1:n_modes
    [~,expModes.q(i)] = max(abs(psiExp_m(:,i)));
    psiExp_m(:,i) = psiExp_m(:,i) / psiExp_m(expModes.q(i),i);
end

expModes.lambdaExp = lambdaExp;
expModes.psiExp = psiExp_m;
expModes.measDOFs = measDOFs;
expModes.lambdaWeights = ones(n_modes,1);
expModes.psiWeights = ones(n_modes,1);

%% Model updating parameter
updatingOpts.formID = 1.0;       % 1: Modal property diff (MAC) ;
                                % 2: Modal property diff (V_mDiff);
updatingOpts.modeMatch = 2;      % 1: Without forced matching;
% 2: With forced matching;
updatingOpts.simModesForExpMatch = 1:n_modes;
updatingOpts.x_lb = -ones(n_alpha,1);
updatingOpts.x_ub =  ones(n_alpha,1);


%% Load local optimal results
load Non_CVX

% Plot alpha value
figHand = figure; set (figHand, 'Position',[250 250 500 250]);
hold on
for j = 1 : 15    
    h = bar(j,local(j),0.15,'k') ;
end

FZ = 10.5;
grid on
ylabel('Relative change \fontname{Times New Roman}{\it\alpha_i}','fontsize',FZ)
xtl = {'{\it\alpha}_1','{\it\alpha}_2','{\it\alpha}_3','{\it\alpha}_4','{\it\alpha}_5','{\it\alpha}_6',...
       '{\it\alpha}_{7}','{\it\alpha}_{8}','{\it\alpha}_{9}','{\it\alpha}_{10}','{\it\alpha}_{11}',...
       '{\it\alpha}_{12}','{\it\alpha}_{13}','{\it\alpha}_{14}','{\it\alpha}_{15}'} ;
h = my_xticklabels(gca,1:15,xtl);
ylim([-1,1]);
hold off


%% Plot relative error

error = abs(local - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
figHand = figure; set (figHand, 'Position',[250 250 500 250]);
hold on
for j = 1 : 15
     h = bar(j,(error(j)),0.15,'k') ; 
end

FZ = 10.5;
grid on
ylabel('Relative error \fontname{Times New Roman}{\ite_i}(%)','fontsize',FZ)
xtl = {'{\itE}_1','{\itE}_2','{\itE}_3','{\itE}_4','{\itE}_5','{\itE}_6',...
    '{\itE}_{t2}','{\itE}_{t3}','{\itE}_{t4}','{\itE}_{t5}','{\itE}_{t6}',...
    '{\itk}_{y1}','{\itk}_{z1}','{\itk}_{y2}','{\itk}_{z2}'} ;
h = my_xticklabels(gca,1:15,xtl);

text(0.2,80,'Obj. Val. = 4.28\times10^{-8}','Fontsize',FZ)
hold off


% Hyperline walk

unmeasDOFs = setdiff( (1 : N)', expModes.measDOFs);
orderIndex = [expModes.measDOFs; unmeasDOFs];
structModel.M0 = structModel.M0(orderIndex, orderIndex);
structModel.K0 = structModel.K0(orderIndex, orderIndex);
for i = 1 : n_alpha
    structModel.K_j{i} = structModel.K_j{i}(orderIndex, orderIndex);
end

fun = @(x)LsqnonlinObjJac(x, structModel, expModes, updatingOpts);


grad = (alpha_act - local)/100;
t = -70:1:215;

fval_step = zeros(length(t),1);

for i = 1:length(t)
    xPoint = local + t(i) * grad;
    f_alpha(i) = norm(feval(fun,xPoint))^2;
    
end

figure
semilogy(t/100, f_alpha,'k','LineWidth',1.5);

xlabel('\theta','FontSize',14);
ylabel('Objective  function \it{f}(\rm{\bf\alpha}(\theta))','FontSize',14);

set(gca,'FontSize',14)

xlim([-0.72 2.20]);
text(0,9e-9,'4.278\times10^{-8}','Fontsize',14)
text(1,1e-29,'7.215\times10^{-29}','Fontsize',14)










