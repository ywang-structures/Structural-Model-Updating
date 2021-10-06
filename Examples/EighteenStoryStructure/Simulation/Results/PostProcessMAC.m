clc;clear;close all

alpha_act = [0.05;  0.05; -0.05; -0.10;  0.10; -0.15;
             0.15;  0.25; -0.10;  0.20;  0.30;  0.25;
            -0.15;  0.05; -0.15;  0.10;  0.20;  0.20;];         
n_alpha = length(alpha_act);
num_star = 100;



%% Plot
filename = 'EighteenStoryStructureSim_form1_JACoff_Levenberg-Marquardt';
load(filename);
fval1 = fval; 
[~,idMin] = min(fval);
time1 = sum(t)+sum(t_waste)
alpha1 = x(1:n_alpha,:);
error(1,:) = abs(alpha1(:,idMin) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
for i = 1:num_star
    error_iter = abs(alpha1(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_error(i,1) = mean(error_iter);
end 

filename = 'EighteenStoryStructureSim_form1_JACon_Levenberg-Marquardt';
load(filename);

fval2 = fval;
[~,idMin] = min(fval);
time2 = sum(t)+sum(t_waste)
alpha2 = x(1:n_alpha,:);
error(2,:) = abs(alpha2(:,idMin) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;

for i = 1:num_star
    error_iter = abs(alpha2(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_error(i,2) = mean(error_iter);
end

% plot;
FZ = 12;         
color = {'k','c','y','w'} ;
k = [-0.10 0.10] ;
figHand = figure; set (figHand, 'Position',[250 250 550 250]);
hold on
for j = 1 : n_alpha
    for i = 1 :2
        h = bar(j + k(i),(error(i,j)),0.15,color{i}) ;
    end
end
grid on
ylabel('Relative error \fontname{Times New Roman}{\ite_i}(%)','fontsize',FZ)
xtl = {'{\itk}_1','{\itk}_2','{\itk}_3','{\itk}_4','{\itk}_5','{\itk}_6',...
       '{\itk}_7','{\itk}_8','{\itk}_9','{\itk}_{10}','{\itk}_{11}','{\itk}_{12}',...
       '{\itk}_{13}','{\itk}_{14}','{\itk}_{15}','{\itk}_{16}','{\itk}_{17}','{\itk}_{18}'};
h = my_xticklabels(gca,1:n_alpha,xtl);
lgd = legend('Numerical Jacobian','Analytical Jacobian');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);
errorAvg = mean(error,2)
box on
ylim([0 6.5])

figHand = figure; 
set (figHand, 'Position',[250 250 550 250]);
plot(1:num_star,mean_error(:,1),'*','Color', [17/255 17/255 17/255])
hold on
box on
plot(1:num_star,mean_error(:,2),'ok')
xlabel('Starting Points','Fontsize',11);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite}_{avg}(%)','fontsize',FZ);
lgd = legend('Numerical Jacobian','Analytical Jacobian');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);
set(gca,'Fontsize',FZ);      