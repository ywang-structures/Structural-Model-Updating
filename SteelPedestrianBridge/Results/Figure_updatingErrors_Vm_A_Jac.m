
% lsqnonlin
         
color = {'k','c','y','w'} ;
k = [-0.10 0.10] ;
figHand = figure; set (figHand, 'Position',[250 250 500 250]);
hold on
for j = 1 : 15
    for i = 1 :2
        h = bar(j + k(i),(error(i,j)),0.15,color{i}) ;
    end
end

FZ = 10.5;
grid on

ylabel('Relative error \fontname{Times New Roman}{\ite_i}(%)','fontsize',FZ)
xtl = {'{\itE}_1','{\itE}_2','{\itE}_3','{\itE}_4','{\itE}_5','{\itE}_6',...
       '{\itE}_{t2}','{\itE}_{t3}','{\itE}_{t4}','{\itE}_{t5}','{\itE}_{t6}',...
       '{\itk}_{y1}','{\itk}_{z1}','{\itk}_{y2}','{\itk}_{z2}'} ;
h = my_xticklabels(gca,1:15,xtl);

legend('Case 2(a)','Case 2(b)')
text(0.2,4.7e-11,'Case 2(a), Avg. Err. = 1.06\times10^{-11}%','Fontsize',FZ)
text(0.2,4.1e-11,'Case 2(b), Avg. Err. = 1.10\times10^{-11}%','Fontsize',FZ)
ylim([0,5e-11])



