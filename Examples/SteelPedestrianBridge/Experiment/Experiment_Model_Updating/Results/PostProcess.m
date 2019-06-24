%%
clc
clear
close all

n_alpha = 17;
FZ = 14;
%% MAC Value formulation

load SteelPedBridg_form1_JAConinterior-point

[~,index_sort] = sort(fval);
figHand = figure;
set (figHand, 'Position',[250 250 450 180]);
semilogy(index_sort,fval(index_sort),'*k');
xlabel('Starting Points #','Fontsize',FZ);
ylabel('Obj. Value','Fontsize',FZ);
set(gca,'fontsize',FZ);
ylim([0,1e10]);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% Plot alpha value
figHand = figure; 
set (figHand, 'Position',[250 250 450 180]);
hold on
for j = 1 : n_alpha
    
    if(j == 15 || j == 17)
        h = bar(j,alpha(j,index_sort(1))/10,0.15,'k');
    else
        h = bar(j,alpha(j,index_sort(1)),0.15,'k');
    end
end

fval(index_sort(1));

grid on
ylabel('{\it\alpha_i^*}','fontsize',FZ,'fontname','Times New Roman')
xtl = {'$\alpha_1$','$\alpha_2$','$\alpha_3$','$\alpha_4$','$\alpha_5$','$\alpha_6$',...
    '$\alpha_7$','$\alpha_8$','$\alpha_9$','$\alpha_{10}$','$\alpha_{11}$','$\alpha_{12}$'...
    ,'$\alpha_{13}$','$\alpha_{14}$','$\frac{\alpha_{15}}{10}$','$\alpha_{16}$','$\frac{\alpha_{17}}{10}$'};
h = my_xticklabelsV2(gca,1 : n_alpha,xtl);
set(gca,'fontsize',FZ-2);
text(0.2,0.6,'Obj. Value = 0.2267','fontsize',FZ)
ylim([-1,1]);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
% 
%% Eigenvector difference formulation

load SteelPedBridg_form2_JAConinterior-point

[~,index_sort] = sort(fval);
figHand = figure;
set (figHand, 'Position',[250 250 450 190]);
semilogy(index_sort,fval(index_sort),'*k');
xlabel('Starting Point #','Fontsize',FZ);
ylabel('Obj. Value','Fontsize',FZ);
set(gca,'fontsize',FZ);
ylim([0,1e4]);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% Plot alpha value
alpha = x;

figHand = figure;
set (figHand, 'Position',[250 250 450 190]);
hold on
for j = 1 : n_alpha
    
    if(j == 15 || j == 17)
        h = bar(j,alpha(j,index_sort(1))/10,0.15,'k');
    else
        h = bar(j,alpha(j,index_sort(1)),0.15,'k');
    end
end

fval(index_sort(1));

grid on
ylabel('{\it\alpha_i^*}','fontsize',FZ,'fontname','Times New Roman')
xtl = {'$\alpha_1$','$\alpha_2$','$\alpha_3$','$\alpha_4$','$\alpha_5$','$\alpha_6$',...
    '$\alpha_7$','$\alpha_8$','$\alpha_9$','$\alpha_{10}$','$\alpha_{11}$','$\alpha_{12}$'...
    ,'$\alpha_{13}$','$\alpha_{14}$','$\frac{\alpha_{15}}{10}$','$\alpha_{16}$','$\frac{\alpha_{17}}{10}$'};
h = my_xticklabelsV2(gca,1 : n_alpha,xtl);
set(gca,'fontsize',FZ);
text(0.2,0.6,'Obj. Value = 2.8474','fontsize',FZ)
ylim([-1,1]);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2)*2;
ax_width = outerpos(3) - ti(1)*0.5;
ax_height = outerpos(4) - ti(2)*2 - ti(4);
ax.Position = [left bottom ax_width ax_height];


%% Dynamic residual

load SteelPedBridg_form3_JACon_trust-region-reflective

[~,index_sort] = sort(fval);
figHand = figure;
set (figHand, 'Position',[250 250 450 190]);
semilogy(index_sort,fval(index_sort),'*k');
xlabel('Starting Point #','Fontsize',FZ);
ylabel('Obj. Value','Fontsize',FZ);
set(gca,'fontsize',FZ);
ylim([3.5e9,4.5e9]);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%% Plot alpha value
alpha = x;
figHand = figure;
set (figHand, 'Position',[250 250 450 190]);
hold on
for j = 1 : n_alpha
    h = bar(j,alpha(j,index_sort(1)),0.15,'k');
end

fval(index_sort(1));


grid on
ylabel('{\it\alpha_i^*}','fontsize',FZ,'fontname','Times New Roman')
xtl = {'$\alpha_1$','$\alpha_2$','$\alpha_3$','$\alpha_4$','$\alpha_5$','$\alpha_6$',...
    '$\alpha_7$','$\alpha_8$','$\alpha_9$','$\alpha_{10}$','$\alpha_{11}$','$\alpha_{12}$'...
    ,'$\alpha_{13}$','$\alpha_{14}$','$\alpha_{15}$','$\alpha_{16}$','$\alpha_{17}$'};

h = my_xticklabelsV2(gca,1 : n_alpha,xtl);
set(gca,'fontsize',FZ);
text(0.2,2.6,'Obj. Value = 3.92\times10^{9}','fontsize',FZ)
ylim([-1,4]);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(4)*1.5;
ax_width = outerpos(3) - ti(1)*0.5;
ax_height = outerpos(4) - ti(4)*2;
ax.Position = [left bottom ax_width ax_height];

