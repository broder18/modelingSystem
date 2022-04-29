function Stat=Statistic(dti,t,q,l,k)
    figure('name',l)
    subplot(2,4,1); hold on; grid minor; %Построение ПВ
    title('ПВ');
    xlabel('{\it X}');
    ylabel('{\it p ( X )}');
    
    r=100;
    if k==1
        r=20000;
    end
    if k==2 
        r=100000;
    end
    if t==1
        hist(dti,r);
    else
        histfit(dti, r, t);
    end
    
    subplot(2,4,2); hold on; grid minor; %Построение ФР
    title('ФР');
    xlabel('{\it X}');
    ylabel('{\it F ( X )}');
    
    [f,x]=ecdf(dti);
    plot(x,f,q);
    
    subplot(2,4,3); hold on; grid minor; %Построение ПВ
    title('ПВ в logY');
    xlabel('{\it X}');
    ylabel('{\it p ( X )}');
    
    r=100;
    if k==1
        r=20000;
    end
    if k==2 
        r=100000;
    end
    if t==1
        hist(dti,r);
    else
        histfit(dti, r, t);
    end
    
    subplot(2,4,4); hold on; grid minor; %Построение ФР
    title('ФР в logX');
    xlabel('{\it X}');
    ylabel('{\it F ( X )}');
    
    [f,x]=ecdf(dti);
    plot(x,f,q);
    
    subplot(2,4,5); hold on; grid minor; %Построение АКФ
    title('Центр. часть АКФ');
    xlabel('{\tau}');
    ylabel('K ( {\tau} )');
    corr=xcorr(dti-mean(dti))/length(dti);
    tau=(-((length(corr)-3)/2)):1:((length(corr)-3)/2);
    corrgr=[corr(1:(((length(corr)-1)/2-1))) max(corr) corr(((length(corr)-1)/2-1):-1:1)];
    plot(tau,corrgr,q);
    xlim([-length(dti)/4 length(dti)/4]);
    
    subplot(2,4,6); hold on; grid minor; %Построение СПМ
    title('СПМ')
    xlabel('{\omega}');
    ylabel('W ( {\omega} )');
    
    S_ft=(abs(fft(corr))/length(corr)).^2;
    plot(1:length(S_ft),S_ft,q);
    
    subplot(2,4,7); hold on; grid minor; %Построение АКФ LOG Y
    title('АКФ в logY');
    xlabel('{\tau}');
    ylabel('K ( {\tau} )');
    
    corr=xcorr(dti-mean(dti))/length(dti);
    tau=(-((length(corr)-3)/2)):1:((length(corr)-3)/2);
    corrgr=[corr(1:(((length(corr)-1)/2-1))) max(corr) corr(((length(corr)-1)/2-1):-1:1)];
    plot(tau,corrgr,q);
    xlim([-length(dti)/4 length(dti)/4]);
    
    subplot(2,4,8); hold on; grid minor; % Построение СПМ LOG Y
    title('СПМ в logY');
    xlabel('{\omega}');
    ylabel('W ( {\omega} )');
    
    S_ft=(abs(fft(corr))/length(corr)).^2;
    plot(1:length(S_ft),S_ft,q);
    
    Stat=1;
end