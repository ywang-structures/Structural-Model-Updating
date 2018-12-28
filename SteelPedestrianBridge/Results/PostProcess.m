clc;clear;close all

alpha_act = [0.05; 0.05; -0.05; -0.10; 0.10; -0.15;
             0.15; -0.05; -0.10; 0.10; -0.20;
             -0.30;0.60;-0.30;0.60;];
         
n_alpha = length(alpha_act);

num_star = 100;

%% MAC Jacon
filename = 'SteelPedBrdg_form1_JACon_Levenberg-Marquardt';
load(filename);
[fval_sort1,index] = sort(fval,'ascend');
x_sort = alpha(:,index);
error(1,:) = abs(x_sort(:,1) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;

for i = 1:num_star
    error_iter = abs(alpha(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_error(i,1) = mean(error_iter);
end

Figure_updatingErrors_MAC_A_Jac



figHand = figure; 
set (figHand, 'Position',[250 250 500 250]);

plot(1:num_star,mean_error(:,1),'*','Color', [17/255 17/255 17/255])
xlabel('Starting Points','Fontsize',FZ);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite}_{avg}(%)','fontsize',FZ);
set(gca,'Fontsize',FZ);
lgd = legend('Case 1(a)');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);


%% Eigval Diff Ana Jac
filename = 'SteelPedBrdg_form2_JACon_Levenberg-Marquardt';
load(filename);
[fval_sort5,index] = sort(fval,'ascend');
x_sort = alpha(:,index);
error(1,:) = abs(x_sort(:,1) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
 
for i = 1:num_star
    error_iter = abs(alpha(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_error(i,1) = mean(error_iter);
end 


filename = 'SteelPedBrdg_form2_JACon_trust-region-reflective';
load(filename);
[fval_sort6,index] = sort(fval,'ascend');
x_sort = alpha(:,index);
error(2,:) = abs(x_sort(:,1) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;


for i = 1:num_star
    error_iter = abs(alpha(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_error(i,2) = mean(error_iter);
end

Figure_updatingErrors_Vm_A_Jac


figHand = figure; 
set (figHand, 'Position',[250 250 500 250]);

plot(1:num_star,mean_error(:,1),'*','Color', [17/255 17/255 17/255])
hold on
plot(1:num_star,mean_error(:,2),'ok')
xlabel('Starting Points','Fontsize',11);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite}_{avg}(%)','fontsize',FZ);
lgd = legend('Case 2(a)','Case 2(b)');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);
set(gca,'Fontsize',FZ);










      