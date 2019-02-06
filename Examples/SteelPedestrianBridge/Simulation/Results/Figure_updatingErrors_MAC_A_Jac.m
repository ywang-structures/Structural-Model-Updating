
% lsqnonlin
         
color = {'k','y','w'} ;
k = [-0.3 -0.10 0.10 0.3] ;
figHand = figure; set (figHand, 'Position',[250 250 500 250]);
hold on
for j = 1 : 15
    for i = 1 :1
        h = bar(j,(error(i,j)),0.15,color{i}) ;
    end
end

FZ = 12;
grid on

ylabel('Relative error \fontname{Times New Roman}{\ite_i}(%)','fontsize',FZ)
xtl = {'{\itE}_1','{\itE}_2','{\itE}_3','{\itE}_4','{\itE}_5','{\itE}_6',...
       '{\itE}_{t2}','{\itE}_{t3}','{\itE}_{t4}','{\itE}_{t5}','{\itE}_{t6}',...
       '{\itk}_{y1}','{\itk}_{z1}','{\itk}_{y2}','{\itk}_{z2}'} ;
h = my_xticklabels(gca,1:15,xtl);

lgd = legend('Case 1(a)');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);
text(0.2,7.5e-3,'Case 1(a), Avg. Err. = 2.00\times10^{-3}%','Fontsize',FZ)
ylim([0,10e-3])



