%%
clear
close all
cd ..
LoadStructure;
cd Results

comCase = 1;
Modeindex = [1 2 3 4 5];

%% Assemble structure parameter:
structModel.M0 = M0;
structModel.K0 = K0;
structModel.K_j = K_j;


%% Experimental modal properties
load ../ExpModeInfo_FEM

numModes = length(Modeindex);
Lambda_e = lambdaExp(Modeindex);
UDOF = setdiff(1:N,measDOFs);

reord = [measDOFs';UDOF'];

K0 = K0(reord,reord);
M0 = M0(reord,reord);

[V,D] = eigs(K0,M0,20,'sm');
[Lambda_nom,srt_idx] = sort(diag(D),'ascend');
V_nom = V(1:length(measDOFs),srt_idx);
V_nom = V_nom(:,Modeindex);
Lambda_nom = Lambda_nom(Modeindex);
f_nom = sqrt(Lambda_nom) / 2/ pi;
MAC_nom = diag(mac(psiExp_m(:,Modeindex),V_nom));

for i = 1 : n_alpha
    K_j{i} = K_j{i}(reord, reord);
end

%% Cross Check
if(comCase == 1)
    load SteelPedBridg_form1_JAConinterior-point
    loc = alpha(:,69);
elseif(comCase == 2)
    load SteelPedBridg_form2_JAConinterior-point
    loc = alpha(:,24);
end


figHand = figure; 
set (figHand, 'Position',[250 250 450 180]);
hold on
for j = 1 : n_alpha
    
    if(j == 15 || j == 17)
        h = bar(j,loc(j)/10,0.15,'k');
    else
        h = bar(j,loc(j),0.15,'k');
    end
end

FZ = 14;
grid on
ylabel('{\it\alpha_i^*}','fontsize',FZ,'fontname','Times New Roman')
xtl = {'$\alpha_1$','$\alpha_2$','$\alpha_3$','$\alpha_4$','$\alpha_5$','$\alpha_6$',...
    '$\alpha_7$','$\alpha_8$','$\alpha_9$','$\alpha_{10}$','$\alpha_{11}$','$\alpha_{12}$'...
    ,'$\alpha_{13}$','$\alpha_{14}$','$\frac{\alpha_{15}}{10}$','$\alpha_{16}$','$\frac{\alpha_{17}}{10}$'};
h = my_xticklabels(gca,1 : n_alpha,xtl);

set(gca,'fontsize',FZ-2);
if(comCase == 1)
    text(0.2,0.6,'Obj. Value = 0.2419','fontsize',FZ)
else
    text(0.2,0.6,'Obj. Value = 2.1448','fontsize',FZ)
end

[~,optIdx] = min(fval);
opt = alpha(:,optIdx);


grad = (opt - loc)/1000;

step_CNT = -10:1:1010;
fval_step = zeros(length(step_CNT),1);

for step = 1:length(step_CNT)
    K_opt = K0;
    xPoint(:,step) = loc + step_CNT(step)*grad;
    for i = 1:n_alpha
        K_opt = K_opt + K_j{i} * xPoint(i,step);
    end
    
    [V,D] = eigs(K_opt,M0,20,'sm');
    [Lambda_opt,srt_idx] = sort(diag(D),'ascend');
    V_opt = V(1:length(measDOFs),srt_idx);
    simModes.psi_u = V(length(measDOFs) + 1 : end,  srt_idx);
    V_opt = V_opt(:,Modeindex);
    simModes.psi_u = simModes.psi_u(:,Modeindex);
    Lambda_opt = Lambda_opt(Modeindex);
    f_opt = sqrt(Lambda_opt) / 2/ pi;
    MAC_opt = diag(mac(psiExp_m(:,Modeindex),V_opt));
    
    
    simModes.Lambda = Lambda_opt;
    simModes.psi_m = V_opt;
    
    expModes.n_meas = length(expModes.measDOFs);
    expModes.n_modes = length(Modeindex);
    
    for i = 1 : expModes.n_modes
        expModes.q(i) = find(abs(expModes.psiExp(:,i)) == 1);
        expModes.psiExp(:,i) = expModes.psiExp(:,i) / expModes.psiExp(expModes.q(i), i);
        simModes.psi_u(:,i) = simModes.psi_u(:,i) / simModes.psi_m(expModes.q(i),i) ;
        simModes.psi_m(:,i) = simModes.psi_m(:,i) / simModes.psi_m(expModes.q(i),i) ;
    end
    
    r = ModelUpdatingObjective(xPoint(:,step), structModel, expModes, simModes, updatingOpts);
    PDiff(step) = norm(r)^2;
    
end


figure
plot(step_CNT/1000, PDiff,'k','LineWidth',1.5);

xlabel('\theta','FontSize',16);
ylabel('Obj.  func. value \it{f}(\rm{\bf\alpha}(\theta))','FontSize',16);

set(gca,'FontSize',16)

%
xlim([-0.10 1.10]);

if(comCase == 1)
    ylim([0.220 0.250]);
    text(0,0.2419,'0.2419','Fontsize',16)
    text(1,0.2267,'0.2267','Fontsize',16)
elseif(comCase == 2)
    ylim([2.1300 2.1550]);
    text(0,2.1448,'2.1448','Fontsize',16)
    text(1,2.1328,'2.1328','Fontsize',16)
    
end

