clc;clear;close all

alpha_act = [0.9;1.2;1.2;0.95;1.2;1.15;
             1.15;1.1;0.9;0.85;1.2;1.15;];             
n_alpha = length(alpha_act);
actual = alpha_act - ones(n_alpha,1);
num_star = 100;
FZ = 12;


%% Modal property difference with MAC value approach, L-M algorithm
filename = 'ConcBuildFrm_form1_JACoff_Levenberg-Marquardt';
load(filename);
[fval_sort1,index] = sort(fval,'ascend');
x_sort = alpha(:,index);
n_alpha = size(alpha,1);
errorMAC(1,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;
for i = 1:num_star
    error_iterA = abs(alpha(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_errorA(i,1) = mean(error_iterA);
end

filename = 'ConcBuildFrm_form1_JACon_Levenberg-Marquardt';
load(filename);
[fval_sort2,index] = sort(fval,'ascend');
x_sort = alpha(:,index);
n_alpha = size(alpha,1);
errorMAC(2,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;
for i = 1:num_star
    error_iterB = abs(alpha(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_errorB(i,1) = mean(error_iterB);
end

figHand = figure; 
set (figHand, 'Position',[250 250 550 250]);
plot(1:num_star,mean_errorA,'*k'); hold on
plot(1:num_star,mean_errorB,'ok')
xlabel('Starting Points','Fontsize',FZ);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite}_{avg} (%)','fontsize',FZ);
lgd = legend('Numerical gradient','Analytical gradient');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);
ylim([0 60])

color = {'k','c','y','w'} ;
k = [-0.10 0.10] ;
figHand = figure; set (figHand, 'Position',[250 250 550 250]);
hold on
for j = 1 : 12
    for i = 1 :2
        h = bar(j + k(i),(errorMAC(i,j)),0.15,color{i}) ;
    end
end
grid on;box on
ylabel('Relative error \fontname{Times New Roman}{\ite_i}(%)','fontsize',FZ)
xtl = {'{\itE}_1','{\itE}_2','{\itE}_3','{\itE}_4','{\itE}_5','{\itE}_6',...
       '{\itE}_7','{\itE}_8','{\itE}_9','{\itE}_{10}','{\itE}_{11}','{\itE}_{12}'};
h = my_xticklabels(gca,1:12,xtl);
lgd = legend('Numerical gradient','Analytical gradient');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);



%% Modal property difference with eigenvector difference approach, TRR algorithm
filename = 'ConcBuildFrm_form2_JACoff_trust-region-reflective';
load(filename);
[fval_sort5,index] = sort(fval,'ascend');
x_sort = alpha(:,index);
errorEig(1,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;
for i = 1:num_star
    error_iterC = abs(alpha(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_errorC(i,1) = mean(error_iterC);
end

filename = 'ConcBuildFrm_form2_JACon_trust-region-reflective';
load(filename);
[fval_sort6,index] = sort(fval,'ascend');
x_sort = alpha(:,index);
errorEig(2,:) = abs(x_sort(:,1) - actual) ./ (ones(n_alpha,1) + actual) * 100;
for i = 1:num_star
    error_iterD = abs(alpha(:,i) - actual) ./ (ones(n_alpha,1) + actual) * 100;
    mean_errorD(i,2) = mean(error_iterD);
end
figHand = figure; 
set (figHand, 'Position',[250 250 550 250]);
plot(1:num_star,mean_errorC,'*k');hold on
plot(1:num_star,mean_errorD,'ok')
xlabel('Starting Points','Fontsize',FZ);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite}_{avg} (%)','fontsize',FZ)
lgd = legend('Numerical gradient','Analytical gradient');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);

color = {'k','c','y','w'} ;
k = [-0.10 0.10] ;
figHand = figure; set (figHand, 'Position',[250 250 550 250]);
hold on
for j = 1 : 12
    for i = 1 :2
        h = bar(j + k(i),(errorEig(i,j)),0.15,color{i}) ;
    end
end
grid on;box on
ylabel('Relative error \fontname{Times New Roman}{\ite_i}(%)','fontsize',FZ)
xtl = {'{\itE}_1','{\itE}_2','{\itE}_3','{\itE}_4','{\itE}_5','{\itE}_6',...
       '{\itE}_7','{\itE}_8','{\itE}_9','{\itE}_{10}','{\itE}_{11}','{\itE}_{12}'};
h = my_xticklabels(gca,1:12,xtl);
lgd = legend('Numerical gradient','Analytical gradient');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);











      