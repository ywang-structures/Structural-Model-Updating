
         
color = {'k','y','w'} ;
k = [-0.3 -0.10 0.10 0.3] ;
figHand = figure; set (figHand, 'Position',[250 250 500 250]);
hold on
for j = 1 : 12
    for i = 1:1
        h = bar(j,(error(i,j)),0.15,color{i}) ;
    end
end


grid on

ylabel('Relative error \fontname{Times New Roman}{\ite_i}(%)','fontsize',FZ)
xtl = {'{\itE}_1','{\itE}_2','{\itE}_3','{\itE}_4','{\itE}_5','{\itE}_6',...
       '{\itE}_7','{\itE}_8','{\itE}_9','{\itE}_{10}','{\itE}_{11}','{\itE}_{12}'};

h = my_xticklabels(gca,1:12,xtl);

lgd = legend('Case 1(a)');
text(0.2,6,'Case 1(a), Avg. Err. = 3.01%','Fontsize',FZ)

set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);




