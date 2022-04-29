function y=Herst1(in,N)
S_tr=in;
% N=input('Показатель степени аппроксимирующего полинома = ');

%Вычисление показателя Херста методом DFA

data_tr = S_tr - mean(S_tr); % Вычисление куммулятивной суммы
data_tr = cumsum(data_tr);
l_max = length(data_tr)/4;
q1 = 1;
q2 = 2;
q3 = 3;
j1 = 1;
figure('name','DFA')

    for b = 1.6:0.5:log(l_max)
        s = [1:floor(exp(b))];% Задание отсчетов окна
        k = length(s);
        Cikl = floor(length(data_tr)/k);
        f1 = 0;
        f2 = 0;
        f3 = 0;
        for i = 1:1:Cikl
            p = polyfit(s-(i+0.5)*k,data_tr(s),N);% Проход по окнам
            z = polyval(p,s-(i+0.5)*k);    
            x = sum((z-data_tr(s)).^2)/k;
            y1 = (x.^(q1/2));
            y2 = (x.^(q2/2));
            y3 = (x.^(q3/2));
            f1 = f1+y1;
            f2 = f2+y2;
            f3 = f3+y3;
            s = s+floor(exp(b));
        end
        
        f1 = ((1/Cikl)*f1).^(1/q1);
        YDFA1(j1) = log(f1/sqrt(k));
        
        plot(log(floor(exp(b))),log(f1/sqrt(k)),'r^-');
        hold on;
  
        f2 = ((1/Cikl)*f2).^(1/q2);
        plot(log(floor(exp(b))),log(f2/sqrt(k)),'ro-'); hold on;
        YDFA2(j1) = log(f2/sqrt(k));

        f3 = ((1/Cikl)*f3).^(1/q3);
        plot(log(floor(exp(b))),log(f3/sqrt(k)),'ro-'); hold on;
        YDFA3(j1) = log(f3/sqrt(k));
        XDFA(j1) = b;
        j1 = j1 + 1;
    end
title('Флуктуационные функции DFA')
xlabel('Величина окна')
ylabel('log(F(q)/d)), d - размер окна')
% Показатель Херста вывод
j=1;
for i = 1:(length(YDFA1))
            XDFA1(j) = XDFA(i);
            YDFA11(j) = YDFA1(i);
            YDFA21(j) = YDFA2(i);
            YDFA31(j) = YDFA3(i);
            j = j + 1;
end

X1 = [ones(length(XDFA1),1) (XDFA1)'];
[b1,bint1] = regress(YDFA11',X1,0.05);
Yy = b1(1)+ b1(2).*XDFA1;
Yn = bint1(1,1)+ bint1(2,1).*XDFA1;
Yv = bint1(1,2)+ bint1(2,2).*XDFA1;
plot(XDFA1,Yy,'r', XDFA1,Yn,'b--', XDFA1,Yv,'g-- '); hold on;
HDFA1=b1(2);
    
X2 = [ones(length(XDFA1),1) (XDFA1)'];
[b2,bint2] = regress(YDFA21',X2,0.05);
Yy = b2(1)+ b2(2).*XDFA1;
Yn = bint2(1,1)+ bint2(2,1).*XDFA1;
Yv = bint2(1,2)+ bint2(2,2).*XDFA1;
plot(XDFA1,Yy,'r', XDFA1,Yn,'b--', XDFA1,Yv,'g-- '); hold on;
HDFA2=b2(2);
    
X3 = [ones(length(XDFA1),1) (XDFA1)'];
[b3,bint3] = regress(YDFA31',X3,0.05);
Yy = b3(1)+ b3(2).*XDFA1;
Yn = bint3(1,1)+ bint3(2,1).*XDFA1;
Yv = bint3(1,2)+ bint3(2,2).*XDFA1;
plot(XDFA1,Yy,'r', XDFA1,Yn,'b--', XDFA1,Yv,'g-- '); hold on;
HDFA3=b3(2);

% --- WTA
n = length(S_tr);
lmax = n/4;
por = 1;
j2 = 1;
figure('name','WTA')

for b = 1.6:0.5:log(lmax)
    s = [1:floor(exp(b))];
    k = length(s);
    Cikl = floor(n/k-por);
    y = zeros(1,Cikl);
    for i = 1:1:Cikl
        y(i) = sum(S_tr(s));
        s = s+floor(exp(b));
    end
    WT = diff(y,por);
    
    WT1 = (abs(WT)).^q1;
    f1 = ((1/Cikl)*sum(WT1)).^(1/q1);
    
    plot(log(floor(exp(b))),log(f1/sqrt(k)),'b^-'); hold on;
    YWTA1(j2) = log(f1/sqrt(k));

    WT2 = (abs(WT)).^q2;
    f2 = ((1/Cikl)*sum(WT2)).^(1/q2);
    plot(log(floor(exp(b))),log(f2/sqrt(k)),'b^-'); hold on;
    YWTA2(j2) = log(f2/sqrt(k));

    WT3 = (abs(WT)).^q3;
    f3 = ((1/Cikl)*sum(WT3)).^(1/q3);
    plot(log(floor(exp(b))),log(f3/sqrt(k)),'b^-'); hold on;
    YWTA3(j2) = log(f3/sqrt(k));

    XWTA(j2) = b;
    j2 = j2 + 1;
end
title('Флуктуационная функция WTA')
xlabel('Величина окна')
ylabel('log(F(q)/d)), d - размер окна')
% WTA вывод
    for i = 1:length(YWTA1)
                XWTA1(j) = XWTA(i);
                YWTA11(j) = YWTA1(i);
                YWTA21(j) = YWTA2(i);
                YWTA31(j) = YWTA3(i);
                j = j + 1;
    end
    
X1 = [ones(length(XWTA1),1) (XWTA1)'];
[b1,bint1] = regress(YWTA11',X1,0.05);
Yy = b1(1)+ b1(2).*XWTA1;
Yn = bint1(1,1)+ bint1(2,1).*XWTA1;
Yv = bint1(1,2)+ bint1(2,2).*XWTA1;
plot(XWTA1,Yy,'r', XWTA1,Yn,'b--', XWTA1,Yv,'g-- '); hold on;
HWTA1=b1(2);
    
X2 = [ones(length(XWTA1),1) (XWTA1)'];
[b2,bint2] = regress(YWTA21',X2,0.05);
Yy = b2(1)+ b2(2).*XWTA1;
Yn = bint2(1,1)+ bint2(2,1).*XWTA1;
Yv = bint2(1,2)+ bint2(2,2).*XWTA1;
plot(XWTA1,Yy,'r',XWTA1,Yn,'b--', XWTA1,Yv,'g-- '); hold on;
HWTA2=b2(2);
    
X3 = [ones(length(XWTA1),1) (XWTA1)'];
[b3,bint3] = regress(YWTA31',X3,0.05);
Yy = b3(1)+ b3(2).*XWTA1;
Yn = bint3(1,1)+ bint3(2,1).*XWTA1;
Yv = bint3(1,2)+ bint3(2,2).*XWTA1;
plot(XWTA1,Yy,'r', XWTA1,Yn,'b--', XWTA1,Yv,'g-- '); hold on;
HWTA3=b3(2);
    
% --- CMA
data1 = S_tr - mean(S_tr);
data1 = cumsum(data1);
l = length(data1);
lmax = l/4;
j3 = 1;
figure('name','CMA')


for b = 1.6:0.5:log(lmax)
    s = [1:floor(exp(b))];
    k = length(s);
    Cikl = floor(l/k);
    y = zeros(1,(2*Cikl)-1);
    for i = 1:1:((2*Cikl)-1)
        y(i) = (sum((data1(s) - mean(data1(s))).^2))/k;
        s = s+floor(exp(b)/2);
    end
    CMA1 = y.^(q1/2);
    f1 = ((1/(2*Cikl))*sum(CMA1)).^(1/q1);
    
    
    plot(log(floor(exp(b))),log(f1/sqrt(k)),'b^-'); hold on;
    YCMA1(j3) = log(f1/sqrt(k));
    
    CMA2 = y.^(q2/2);
    f2 = ((1/(2*Cikl))*sum(CMA2)).^(1/q2);
    plot(log(floor(exp(b))),log(f2/sqrt(k)),'b^-'); hold on;
    YCMA2(j3) = log(f2/sqrt(k));
        
    CMA3 = y.^(q3/2);
    f3 = ((1/(2*Cikl))*sum(CMA3)).^(1/q3);
    plot(log(floor(exp(b))),log(f3/sqrt(k)),'b^-'); hold on;
    YCMA3(j3) = log(f3/sqrt(k));
        
    XCMA(j3) = b;
    j3 = j3 + 1;
end
xlabel('Величина окна')
ylabel('log(F(q)/d)), d - размер окна')
title('Флуктуационные функции CMA')
%CMA вывод
    for i = 1:length(YCMA1)
                XCMA1(j) = XCMA(i);
                YCMA11(j) = YCMA1(i);
                YCMA21(j) = YCMA2(i);
                YCMA31(j) = YCMA3(i);
                j = j + 1;
    end
X1 = [ones(length(XCMA1),1) (XCMA1)'];
[b1,bint1] = regress(YCMA11',X1,0.05);
Yy = b1(1)+ b1(2).*XCMA1;
Yn = bint1(1,1)+ bint1(2,1).*XCMA1;
Yv = bint1(1,2)+ bint1(2,2).*XCMA1;
plot(XCMA1,Yy,'r', XCMA1,Yn,'b--', XCMA1,Yv,'g-- '); hold on;
HCMA1=b1(2);

X2 = [ones(length(XCMA1),1) (XCMA1)'];
[b2,bint2] = regress(YCMA21',X2,0.05);
Yy = b2(1)+ b2(2).*XCMA1;
Yn = bint2(1,1)+ bint2(2,1).*XCMA1;
Yv = bint2(1,2)+ bint2(2,2).*XCMA1;
plot(XCMA1,Yy,'r', XCMA1,Yn,'b--', XCMA1,Yv,'g-- '); hold on;
HCMA2=b2(2);
    
X3 = [ones(length(XCMA1),1) (XCMA1)'];
[b3,bint3] = regress(YCMA31',X3,0.05);
Yy = b3(1)+ b3(2).*XCMA1;
Yn = bint3(1,1)+ bint3(2,1).*XCMA1;
Yv = bint3(1,2)+ bint3(2,2).*XCMA1;
plot(XCMA1,Yy,'r', XCMA1,Yn,'b--', XCMA1,Yv,'g-- '); hold on;
HCMA3=b3(2);

Vyh=[HDFA1 HWTA1 HCMA1; HDFA2 HWTA2 HCMA2; HDFA3 HWTA3 HCMA3];
y=Vyh;
end