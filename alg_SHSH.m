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
disp(['Гениальная загадка от Жака Фреско, на компиляцию даётся: ' num2str(count) '0 ' ' минут']); %отображение значений счетчика, чтобы понимать, что все "ОК"
z=S_tek_1;

end