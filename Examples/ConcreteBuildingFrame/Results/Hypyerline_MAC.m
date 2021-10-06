% Updating for NEES: frequency domain
clear
close all
idxLocal = 97
alpha_act = [-0.1;  0.2;   0.2; -0.05; 0.2;  0.15;
             0.15; 0.10; -0.10; -0.15; 0.20; 0.15;];
n_alpha = length(alpha_act);
n_modes = 3;
modeIndex = 1:n_modes;
num_star = 100;
FZ = 12;



%% Modal property difference with MAC value approach, L-M algorithm
filename = 'ConcBuildFrm_form1_JACon_Levenberg-Marquardt';
load(filename);
[fval_sortB,index] = sort(fval,'ascend');
x_sort = alpha(:,index);
n_alpha = size(alpha,1);
errorMAC(2,:) = abs(x_sort(:,1) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
for i = 1:num_star
    error_iterB = abs(alpha(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_errorB(i,1) = mean(error_iterB);
end

filename = 'ConcBuildFrm_form1_JACoff_Levenberg-Marquardt';
load(filename);
[fval_sortA,index] = sort(fval,'ascend');
x_sort = alpha(:,index);
grad_sortA = gradient(:,index);
n_alpha = size(alpha,1);
errorMAC(1,:) = abs(x_sort(:,1) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
for i = 1:num_star
    error_iterA = abs(alpha(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_errorA(i,1) = mean(error_iterA);
end

figHand = figure; 
set (figHand, 'Position',[250 250 500 250]);
plot(1:num_star,mean_errorA,'*k'); hold on
plot(1:num_star,mean_errorB,'ok')
plot(idxLocal,mean_errorA(idxLocal),'*r');
xlabel('Starting Points','Fontsize',FZ);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite}_{avg} (%)','fontsize',FZ);
lgd = legend('Numerical gradient','Analytical gradient');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);
text(idxLocal-4,mean_errorA(idxLocal),'A','Fontsize',12,'Color','r')
ylim([0 60])



%% Assemble structure paramter
M0 = structModel.M0;
K0 = structModel.K0;
K_j = structModel.K_j;
K_act = K0;
for i = 1:n_alpha
    K_act = K_act + K_j{i} * alpha_act(i);
end
N = size(M0,2);
measDOFs = expModes.measDOFs';



%% Simulate "experimental data"
[psiExp,lambdaExp] = eigs(K_act,M0,20,'sm');
[lambdaExp,dummyInd] = sort((diag(lambdaExp)),'ascend') ;
lambdaExp = lambdaExp(modeIndex);
psiExp = psiExp(:,dummyInd(modeIndex));
psiExp_m = psiExp(measDOFs,:);


psiExp_m = expModes.psiExp;
lambdaExp = expModes.lambdaExp;
for i = 1:n_modes
    [~,expModes.q(i)] = max(abs(psiExp_m(:,i)));
    psiExp_m(:,i) = psiExp_m(:,i) / psiExp_m(expModes.q(i),i);
end



%% Load local optimal results
n_alpha = size(alpha,1);
for i = 1:length(alpha)
    error_iterA = abs(alpha(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_errorA(i,1) = mean(error_iterA);
end

local = alpha(:,idxLocal);
disp(['The objective function value of the local solution is ' num2str(fval(idxLocal))])
disp(['The ave. relative error of the local solution is ' num2str(mean_errorA(idxLocal)) '%'])




%% Hyperline walk
unmeasDOFs = setdiff( (1 : N)', expModes.measDOFs);
orderIndex = [measDOFs; unmeasDOFs];
structModel.M0 = structModel.M0(orderIndex, orderIndex);
structModel.K0 = structModel.K0(orderIndex, orderIndex);
for i = 1 : n_alpha
    structModel.K_j{i} = structModel.K_j{i}(orderIndex, orderIndex);
end
fun = @(x)OptmzObjJac(x, structModel, expModes, updatingOpts, 'lsqnonlin');
grad = (alpha_act - local)/100;
t = -70:1:215;
fval_step = zeros(length(t),1);

for i = 1:length(t)
    xPoint = local + t(i) * grad;
    f_alpha(i) = norm(feval(fun,xPoint))^2;    
end

figHand = figure; set (figHand, 'Position',[250 250 550 300]);
semilogy(t/100, f_alpha,'k','LineWidth',1.5);
xlabel('\theta','FontSize',11);
ylabel('Objective  function \it{f}(\rm{\bf\alpha}(\theta))','FontSize',11);
set(gca,'FontSize',11)
xlim([-0.72 2.20]);
text(-0.034,10^-10,'A','Fontsize',12,'Color','k')
text(0.65,10^-19,'Global optimum','Fontsize',12,'Color','k')