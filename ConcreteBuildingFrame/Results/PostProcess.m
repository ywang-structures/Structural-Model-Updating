clc;clear;close all

Opt_stiff = [0.9;1.2;1.2;0.95;1.2;1.15;
             1.15;1.1;0.9;0.85;1.2;1.15;];
             
n_alpha = length(Opt_stiff);
actual = Opt_stiff - ones(n_alpha,1);

num_star = 100;

FZ = 12;
%% MAC JACoff  LM
filename = 'ConcBuildFrm_form1_JACoff_Levenberg-Marquardt';
load(filename);
[fval_sort1,index] = sort(fval,'ascend');
x_sort = alpha(:,index);
n_alpha = size(alpha,1);
error(1,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;

for i = 1:num_star
    error_iter = abs(alpha(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_error(i,1) = mean(error_iter);
end

figHand = figure; 
set (figHand, 'Position',[250 250 500 250]);

plot(1:num_star,mean_error,'*','Color', [17/255 17/255 17/255])
xlabel('Starting Points','Fontsize',FZ);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite}_{avg} (%)','fontsize',FZ)
lgd = legend('Case 1(a)');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);

Figure_updatingErrors_MAC_N_Jac

%% MAC JACon
filename = 'ConcBuildFrm_form1_JACon_Levenberg-Marquardt';
load(filename);
[fval_sort2,index] = sort(fval,'ascend');
x_sort = alpha(:,index);
n_alpha = size(alpha,1);
error(1,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;
% 
for i = 1:num_star
    error_iter = abs(alpha(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_error(i,1) = mean(error_iter);
end
% 
figHand = figure; 
set (figHand, 'Position',[250 250 500 250]);

plot(1:num_star,mean_error,'*','Color', [17/255 17/255 17/255])
xlabel('Starting Points','Fontsize',11);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite}_{avg} (%)','fontsize',FZ)
lgd = legend('Case 1(a)');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);

Figure_updatingErrors_MAC_A_Jac


%% Eigvec Jacoff

filename = 'ConcBuildFrm_form2_JACoff_Levenberg-Marquardt';
load(filename);
[fval_sort3,index] = sort(fval,'ascend');
x_sort = alpha(:,index);
n_alpha = size(alpha,1);
error(1,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;

for i = 1:num_star
    error_iter = abs(alpha(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_error(i,1) = mean(error_iter);
end

filename = 'ConcBuildFrm_form2_JACoff_trust-region-reflective';
load(filename);
[fval_sort4,index] = sort(fval,'ascend');
x_sort = alpha(:,index);
error(2,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;

for i = 1:num_star
    error_iter = abs(alpha(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_error(i,2) = mean(error_iter);
end


Figure_updatingErrors_Vm_N_Jac


figHand = figure; 
set (figHand, 'Position',[250 250 500 250]);

plot(1:num_star,mean_error(:,1),'*','Color', [17/255 17/255 17/255])
hold on
plot(1:num_star,mean_error(:,2),'ok')
xlabel('Starting Points','Fontsize',11);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite}_{avg} (%)','fontsize',11)
lgd = legend('Case 2(a)','Case 2(b)');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);

%% Eigvec Jacon

filename = 'ConcBuildFrm_form2_JACon_Levenberg-Marquardt';
load(filename);
[fval_sort5,index] = sort(fval,'ascend');
x_sort = alpha(:,index);
error(1,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;

for i = 1:num_star
    error_iter = abs(alpha(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_error(i,1) = mean(error_iter);
end

filename = 'ConcBuildFrm_form2_JACon_trust-region-reflective';
load(filename);
[fval_sort6,index] = sort(fval,'ascend');
x_sort = alpha(:,index);

error(2,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;

for i = 1:num_star
    error_iter = abs(alpha(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_error(i,2) = mean(error_iter);
end

Figure_updatingErrors_Vm_A_Jac


figHand = figure; 
set (figHand, 'Position',[250 250 500 250]);

plot(1:num_star,mean_error(:,1),'*','Color', [17/255 17/255 17/255])
hold on
plot(1:num_star,mean_error(:,2),'ok')
xlabel('Starting Points','Fontsize',11);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite}_{avg} (%)','fontsize',11)
lgd = legend('Case 2(a)','Case 2(b)');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);














      