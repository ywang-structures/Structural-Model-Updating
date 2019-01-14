      
color = {'k','y','w'} ;
k = [-0.3 -0.10 0.10 0.3] ;
figHand = figure; set (figHand, 'Position',[250 250 450 180]);
hold on
for j = 1 : n_alpha
    for i = 1 :1
        h = bar(j,(error(i,j)),0.15,color{i}) ;
    end
end

FZ = 14;
grid on

ylabel('{\it\alpha_i^*}','fontsize',FZ,'fontname','Times New Roman')
xtl = {'{\itk}_1','{\itk}_2','{\itk}_3','{\itk}_{4}'} ;
h = my_xticklabels(gca,1:n_alpha,xtl);

objV = fval_sort1(1);
power = ceil(log10(objV)-1);
base = 10^(log10(objV)-power);

set(gca,'fontsize',FZ);
text(0.2,7.5e-3,sprintf('Obj. Val. = %0.2f\\times10^{%d}',base,power),'Fontsize',FZ)
ylim([-1,1])



