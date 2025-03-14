function [] = plotting_p(a,b,c) % pearson correlation

    R= corrcoef(a);
    s1 = plot(a(:,1), a(:,2), 'r+');
    set(s1, 'MarkerSize', 8, 'LineWidth', 2);
    s1 = ancestor(s1, 'axes');
    s1.YAxis.Exponent = 0;
    s1.XAxis.Exponent = 0;
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
end  