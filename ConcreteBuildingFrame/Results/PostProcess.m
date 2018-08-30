clc;clear;close all

Opt_stiff = [0.9;1.2;1.2;0.95;1.2;1.15;
             1.15;1.1;0.9;0.85;1.2;1.15;];
             
n_alpha = length(Opt_stiff);
actual = Opt_stiff - ones(n_alpha,1);

num_start = 100;


%% MAC JACoff  LM
filename = 'ConcBuildFrm_form1_JACoff_LM';
load(filename);
[fval_sort1,index] = sort(fval,'ascend');
x_sort = x(:,index);
n_alpha = size(x,1);
error(1,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;

for i = 1:num_start
    error_iter = abs(x(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_error(i,1) = mean(error_iter);
end

figHand = figure; 
set (figHand, 'Position',[250 250 500 250]);

plot(1:num_start,mean_error,'ok')
xlabel('Starting Points','Fontsize',11);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite}_{avg} (%)','fontsize',11)
legend('Case 1(a)')

Figure_updatingErrors_MAC_N_Jac

%% MAC JACon
filename = 'ConcBuildFrm_form1_JACon_LM';
load(filename);
[fval_sort2,index] = sort(fval,'ascend');
x_sort = x(:,index);
n_alpha = size(x,1);
error(1,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;
% 
for i = 1:num_start
    error_iter = abs(x(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_error(i,1) = mean(error_iter);
end
% 
figHand = figure; 
set (figHand, 'Position',[250 250 500 250]);

plot(1:num_start,mean_error,'ok')
xlabel('Starting Points','Fontsize',11);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite}_{avg} (%)','fontsize',11)
legend('Case 1(a)')

Figure_updatingErrors_MAC_A_Jac





%% Eigvec Jacoff

filename = 'ConcBuildFrm_form2_JACoff_LM';
load(filename);
[fval_sort3,index] = sort(fval,'ascend');
x_sort = x(:,index);
n_alpha = size(x,1);
error(1,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;

for i = 1:num_start
    error_iter = abs(x(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_error(i,1) = mean(error_iter);
end

filename = 'ConcBuildFrm_form2_JACoff_TRR';
load(filename);
[fval_sort4,index] = sort(fval,'ascend');
x_sort = x(:,index);
n_alpha = size(x,1);
error(2,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;

for i = 1:num_start
    error_iter = abs(x(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_error(i,2) = mean(error_iter);
end


Figure_updatingErrors_Vm_N_Jac


figHand = figure; 
set (figHand, 'Position',[250 250 500 250]);

plot(1:num_start,mean_error(:,1),'*','Color', [17/255 17/255 17/255])
hold on
plot(1:num_start,mean_error(:,2),'ok')
xlabel('Starting Points','Fontsize',11);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite}_{avg} (%)','fontsize',11)
legend('Case 2(a)','Case 2(b)')


%% Eigvec Jacon

filename = 'ConcBuildFrm_form2_JACon_LM';
load(filename);
[fval_sort5,index] = sort(fval,'ascend');
x_sort = x(:,index);
n_alpha = size(x,1);
error(1,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;

for i = 1:num_start
    error_iter = abs(x(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_error(i,1) = mean(error_iter);
end



filename = 'ConcBuildFrm_form2_JACon_TRR';
load(filename);
[fval_sort6,index] = sort(fval,'ascend');
x_sort = x(:,index);
n_alpha = size(x,1);
error(2,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;

for i = 1:num_start
    error_iter = abs(x(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_error(i,2) = mean(error_iter);
end

Figure_updatingErrors_Vm_A_Jac


figHand = figure; 
set (figHand, 'Position',[250 250 500 250]);

plot(1:num_start,mean_error(:,1),'*','Color', [17/255 17/255 17/255])
hold on
plot(1:num_start,mean_error(:,2),'ok')
xlabel('Starting Points','Fontsize',11);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite}_{avg} (%)','fontsize',11)
legend('Case 2(a)','Case 2(b)')


















      