function r=build_cdf(dti,q) %Построение ФР
    [f,x]=ecdf(dti);
    plot(x,f,q);
    xlim([0 length(dti)/3]);
    r=1;
end