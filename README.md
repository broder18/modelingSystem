# Modeling of an infocommunication queuing system based on queue theory

## Функция алгоритма ***Шрайбера-Шмитца***
```
function z=alg_SHSH(S_tek, H, n,N) % Алгоритм Шрайбера-Шмитца
H_tr=H; % Показатель Херста
gamma_tr=2-2*H_tr;
beta_tr=1-gamma_tr;
k_test=1; % начальное условие Стат-теста
count=0; %счетчик итераций
Fd=1; % Частота дискретизации
w=0.00001:2*pi*Fd:((n-1)*Fd*2*pi+0.00001); %сетка частот
S_tek_1=S_tek;
while k_test==1
    count = count+1; % счетчик (инкремент) итераций
    H_tek=Herst(S_tek_1-mean(S_tek_1),N); %Текущий спектр
    H_tek=H_tek(2);
    delta_beta=beta_tr-(1-(2-2*H_tek));
    K2=w.^(-delta_beta/2);
    
    S_tek_fft=fft(S_tek_1-mean(S_tek_1));
    S_tek_fft_vyh=S_tek_fft.*K2;
    S_tek_vyh =real(ifft(S_tek_fft_vyh));
    
    S_tek_vyh=(S_tek_vyh-mean(S_tek_vyh))*std(S_tek)/std(S_tek_vyh)+mean(S_tek);%нормировка
    
    k_test=kstest2(S_tek, S_tek_vyh); %Стат. тест Колмогорова-Смирнова
    %Поранговая замена отсчетов
    [zn1, ~]=sort(S_tek);
    [~, ind2]=sort(S_tek_vyh);
    S_tek_1(ind2)=zn1;
    
    if count==500 %больше 500 итераций не делаем
         return
    end
end
disp(num2str(count)); %отображение значений счетчика, чтобы понимать, что все "ОК"
z=S_tek_1;
end
```
---
## Функция для расчета вероятностей отказа в СМО с отказом
```
function [Potk] = func_CMO_otk(vi, ti, c, U, t) % СМО с отказами
P = zeros(1, length(c));
for i = 1 : 1 : length(c)
    To = ti(1) + vi(1)/c(i); %Время обслуживания 1-го запроса
    for j = 2 : 1 : length(ti) % проверка доступа для текущего запроса
        if ti(j) < To % если время прихода j-го запроса меньше времени обслуживания текущего запроса
            P(i) = P(i) + 1; % то увеличиваем число отказов на один 
        else
            To = ti(j) + vi(j)/c(i); % принимаем запрос и записываем новое время, когда канал будет свободен 
        end
    end
end
 
Potk = P./length(ti); % Матрица вероятностей отказа для различных пропускных способностей канала
 
subplot(2,1,1);
plot(c, Potk, '--', 'LineWidth', 1); xlabel('C'); ylabel('P_{ОТКАЗ} ( C )'); grid on; hold on
subplot(2,1,2);
plot(U, Potk, '--', 'LineWidth', 1); xlabel('U'); ylabel('P_{ОТКАЗ} ( U )'); grid on; hold on 
 
end
```
---
## Функция для расчета среднего времени ожидания и средней длины очереди в СМО с ожиданием
```
function [T_wait_mean L_wait_mean] = func_CMO_wait(vi, ti, c, U, t) % СМО с ожиданиями
L_wait = zeros(1, length(c)); % Средняя длина очереди
T_wait = zeros(1, length(c)); % Среднее время ожидания
 
for j = 1 : 1 : length(c)
    To = ti(1) + vi(1)/c(1); % Время обслуживания 1-го запроса 
    for i = 2 : 1 : length(ti)
        if ti(i) < To % условие попадания в очередь
            T_wait(j) = T_wait(j) + To - ti(i); %время ожидания
            for k = i : 1 : length(ti) % подсчет очереди
                if ti(k) < To
                    L_wait(j) = L_wait(j) + 1;
                else
                    break
                end
            end
            To = To + vi(i)/c(j); % Новое время, когда канал будет свободен, если сформировалась очередь из запросов
        else
            To = ti(i) + vi(i)/c(j); % Новое время, когда канал будет свободен, если в очереди запросов нет
        end
    end
end
 
L_wait_mean = L_wait/length(ti); % Средняя длина очереди - итог
T_wait_mean = T_wait/length(ti); % Среднее время ожидания - итог
 
subplot(2,2,1);
semilogy(c, T_wait_mean, '--', 'LineWidth', 2); xlabel('C'); ylabel('T_{WAIT} ( C )'); grid on; hold on;
subplot(2,2,2);
semilogy(U, T_wait_mean, '--', 'LineWidth', 2); xlabel('U'); ylabel('T_{WAIT} ( U )'); grid on; hold on;
subplot(2,2,3);
semilogy(c, L_wait_mean, '--', 'LineWidth', 2); xlabel('C'); ylabel('L_{WAIT} ( C )'); grid on; hold on;
subplot(2,2,4);
semilogy(U, L_wait_mean, '--', 'LineWidth', 2); xlabel('U'); ylabel('L_{WAIT} ( U )'); grid on; hold on;
 
end
```
## Функция для расчёта плотности вероятности, функции распределения, автокорреляционной функции и спектральной плотности мощности выборки данных
```
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
    semilogy(tau,corrgr,q);
    xlim([-length(dti)/4 length(dti)/4]);
    
    subplot(2,4,8); hold on; grid minor; % Построение СПМ LOG Y
    title('СПМ в logY');
    xlabel('{\omega}');
    ylabel('W ( {\omega} )');
    
    S_ft=(abs(fft(corr))/length(corr)).^2;
    semilogy(1:length(S_ft),S_ft,q);
    
    Stat=1;
end
```
