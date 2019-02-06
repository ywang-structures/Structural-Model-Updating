clc;clear;close all

alpha_act = [-0.1];
         
n_alpha = length(alpha_act);

num_star = 1000;

FZ = 14;

%% MAC Jacon
filename = 'Sim_form1_JACon_interior-point';
load(filename);
[fval_sort1,index] = sort(fval,'ascend');
x = alpha;
x_sort = alpha(:,index);
error(1,:) = abs(x_sort(:,1) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;

for i = 1:num_star
    error_iter = abs(x(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_error(i,1) = mean(error_iter);
end

figHand = figure; 
set (figHand, 'Position',[250 250 450 180]);

plot(1:num_star,mean_error(:,1),'*','Color', [17/255 17/255 17/255])
xlabel('Starting Points','Fontsize',11);
ylabel({'Relative error' 'of \alpha (%)'},'fontsize',FZ);
set(gca,'FontSize',FZ);




%% Eigval Diff Ana Jac
filename = 'Sim_form2_JACon_interior-point';
load(filename);
[fval_sort1,index] = sort(fval,'ascend');
x = alpha;
x_sort = alpha(:,index);
error(1,:) = abs(x_sort(:,1) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;

for i = 1:num_star
    error_iter = abs(x(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_error(i,1) = mean(error_iter);
end

figHand = figure; 
set (figHand, 'Position',[250 250 450 180]);

plot(1:num_star,mean_error(:,1),'*','Color', [17/255 17/255 17/255])
xlabel('Starting Points','Fontsize',11);
ylabel({'Relative error' 'of \alpha (%)'},'fontsize',FZ);
set(gca,'FontSize',FZ);


%% Modal Dynamic Residual Ana Jac
filename = 'Sim_form3_JACon_interior-point';
load(filename);
x = alpha(1:n_alpha,:);
[fval_sort7,index] = sort(fval,'ascend');
x_sort = x(:,index);
n_alpha = size(x,1);
error(1,:) = abs(x_sort(:,1) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
 
for i = 1:num_star
    error_iter = abs(x(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_error(i,1) = mean(error_iter);
end 

figHand = figure; 
set (figHand, 'Position',[250 250 450 180]);

plot(1:num_star,mean_error(1:num_star,1),'*','Color', [17/255 17/255 17/255])
xlabel('Starting Points','Fontsize',11);
ylabel({'Relative error' 'of \alpha (%)'},'fontsize',FZ);
set(gca,'FontSize',FZ);












      