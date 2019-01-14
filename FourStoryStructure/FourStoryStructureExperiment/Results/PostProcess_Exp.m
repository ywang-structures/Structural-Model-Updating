clc;clear;close all
   
n_alpha = 4;

num_star = 1000;

FZ = 14;

%% MAC Jacon
filename = 'Exp_form1_Jacon_interior-point';
load(filename);
[fval_sort1,index] = sort(fval,'ascend');
x = alpha(1:n_alpha,:);
x_sort = alpha(1 : n_alpha,index);
error(1,:) = x_sort(:,1);

Figure_updatingErrors_MP_MAC


figHand = figure; 
set (figHand, 'Position',[250 250 450 180]);

plot(1:num_star,fval,'*','Color', [17/255 17/255 17/255])

hx = xlabel('Starting Points #');
hy = ylabel('Obj. Value');
set(hx, 'FontSize', FZ);
set(hy, 'FontSize', FZ);
set(gca,'Fontsize',FZ);


%% Eigval Diff Ana Jac
filename = 'Exp_form2_Jacon_interior-point';
load(filename);
[fval_sort2,index] = sort(fval,'ascend');
x = alpha(1:n_alpha,:);
x_sort = alpha(1 : n_alpha,index);
error(1,:) = x_sort(:,1);

Figure_updatingErrors_MP_EigDiff

figHand = figure; 
set (figHand, 'Position',[250 250 450 180]);

plot(1:num_star,fval(1:num_star),'*','Color', [17/255 17/255 17/255])
ylim([2.055e-3 2.06e-3])

hx = xlabel('Starting Points #','Fontsize',FZ);
hy = ylabel('Obj. Value','fontsize',FZ);
set(hx, 'FontSize', FZ);
set(hy, 'FontSize', FZ);
set(gca,'Fontsize',FZ);

%% Modal Dynamic Residual
filename = 'Exp_form3_Jacon_interior-point';
load(filename);
[fval_sort3,index] = sort(fval,'ascend');
x = alpha(1:n_alpha,:);
x_sort = alpha(1 : n_alpha,index);
error(1,:) = x_sort(:,1);

Figure_updatingErrors_MDR

figHand = figure; 
set (figHand, 'Position',[250 250 450 180]);

plot(1:num_star,fval(1:num_star),'*','Color', [17/255 17/255 17/255])

hx = xlabel('Starting Points #','Fontsize',FZ);
hy = ylabel('Obj. Value','Fontsize',FZ);
set(hx, 'FontSize', FZ);
set(hy, 'FontSize', FZ);
set(gca,'Fontsize',FZ);






      