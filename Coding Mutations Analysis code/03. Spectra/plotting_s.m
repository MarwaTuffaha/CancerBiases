function [] = plotting_s(a,b,c)

    R= corr(a,'Type','Spearman');
    
    [~,p1] = sort(a(:,1),'ascend');
    r1 = 1:length(a(:,1));
    r1(p1) = r1;

    [~,p2] = sort(a(:,2),'ascend');
    r2 = 1:length(a(:,2));
    r2(p2) = r2;


    s1 = plot(r1, r2, 'b+');
    set(s1, 'MarkerSize', 8, 'LineWidth', 2);
    s1 = ancestor(s1, 'axes');
    s1.YAxis.Exponent = 0;
    
    %%% regression line
    hold on
    l = lsline ;
    set(l,'LineWidth', 2)
    %%% axis display 
    xlabel(b, 'FontSize', 10)
    ylabel(c, 'FontSize', 10)
    set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
    %t=text(min(a(:,1)),max(a(:,2)),{strcat('CorCoef = ',string(R(1,2)))});
    %t.FontSize=15;
    title(strcat("CorCoef=",string(R(1,2))))
    pbaspect([1 1 1])
end   
