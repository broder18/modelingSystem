function y=Herst(in,N)
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
%         figure(2)
%         title('Флуктуационные функции DFA')
%         plot(log(floor(exp(b))),log(f1/sqrt(k)),'r^-');
%         hold on;
  
        f2 = ((1/Cikl)*f2).^(1/q2);
%         plot(log(floor(exp(b))),log(f2/sqrt(k)),'ro-');
        YDFA2(j1) = log(f2/sqrt(k));

        f3 = ((1/Cikl)*f3).^(1/q3);
%         plot(log(floor(exp(b))),log(f3/sqrt(k)),'ro-');
        YDFA3(j1) = log(f3/sqrt(k));
        XDFA(j1) = b;
        j1 = j1 + 1;
    end
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
% figure(2)
% hold on;
% Yy = b1(1)+ b1(2).*XDFA1;
% Yn = bint1(1,1)+ bint1(2,1).*XDFA1;
% Yv = bint1(1,2)+ bint1(2,2).*XDFA1;
% plot(XDFA1,Yy,'r', XDFA1,Yn,'b--', XDFA1,Yv,'g-- ');
HDFA1=b1(2);
    
X2 = [ones(length(XDFA1),1) (XDFA1)'];
[b2,bint2] = regress(YDFA21',X2,0.05);
% Yy = b2(1)+ b2(2).*XDFA1;
% Yn = bint2(1,1)+ bint2(2,1).*XDFA1;
% Yv = bint2(1,2)+ bint2(2,2).*XDFA1;
% plot(XDFA1,Yy,'r', XDFA1,Yn,'b--', XDFA1,Yv,'g-- ');
HDFA2=b2(2);
    
X3 = [ones(length(XDFA1),1) (XDFA1)'];
[b3,bint3] = regress(YDFA31',X3,0.05);
% Yn = bint2(1,1)+ bint2(2,1).*XDFA1;
% Yv = bint2(1,2)+ bint2(2,2).*XDFA1;
% plot(XDFA1,Yy,'r', XDFA1,Yn,'b--', XDFA1,Yv,'g-- ');
HDFA3=b3(2);


Vyh=[HDFA1 HDFA2 HDFA3];
y=Vyh;
end