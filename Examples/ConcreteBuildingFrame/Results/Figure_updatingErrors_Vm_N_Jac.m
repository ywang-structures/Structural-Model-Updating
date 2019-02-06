
% lsqnonlin
         
color = {'k','c','y','w'} ;
k = [-0.10 0.10] ;
figHand = figure; set (figHand, 'Position',[250 250 500 250]);
hold on
for j = 1 : 12
    for i = 1 :2
        h = bar(j + k(i),(error(i,j)),0.15,color{i}) ;
    end
end
set(gca,'yscale','log')

grid on

ylabel('Relative error \fontname{Times New Roman}{\ite_i}(%)','fontsize',FZ)
xtl = {'{\itE}_1','{\itE}_2','{\itE}_3','{\itE}_4','{\itE}_5','{\itE}_6',...
       '{\itE}_7','{\itE}_8','{\itE}_9','{\itE}_{10}','{\itE}_{11}','{\itE}_{12}'};

h = my_xticklabels(gca,1:12,xtl);


lgd = legend('Case 2(a)','Case 2(b)');

set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);

text(0.2,200,'Case 2(a), Avg. Err. = 1.75\times10^{-7}%','Fontsize',FZ)
text(0.2,8,'Case 2(b), Avg. Err. = 0.0165%','Fontsize',FZ)
ylim([0,200])