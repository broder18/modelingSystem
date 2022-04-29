function [T_wait_mean, L_wait_mean] = CMO_ozhid(vi, ti, c, U, t) % СМО с ожиданиями
figure(5);
% Подготовка массивов
L_wait = zeros(1, length(c));
T_wait = zeros(1, length(c)); 
T_obsl = zeros(1, length(c));
T_being = zeros(1,length(c));
parfor j = 1 : 1 : length(c)
    To = ti(1) + vi(1)/c(1); % время обслуживания 1-го запроса
    for i = 2 : 1 : length(ti)
        if ti(i) < To % если попал в очередь
            T_wait(j) = T_wait(j) + To - ti(i); % ждет столько то
            for k = i : 1 : length(ti) % подсчет очереди
                if ti(k) < To
                    L_wait(j) = L_wait(j) + 1;
                end
            end
            To = To + vi(i)/c(j); % Новое время, когда канал будет свободен, если есть очередь 
        else
            To = ti(i) + vi(i)/c(j); % Новое время, когда канал будет свободен, если в очереди никого не было
        end
    T_obsl(j)=T_obsl(j)+vi(i)/c(j);
    T_being(j)=T_being(j)+T_wait(j)+vi(i)/c(j);
    end
end
L_wait_mean = L_wait/length(ti); % Cредняя длина очереди
T_wait_mean = T_wait/length(ti); % Среднее время ожидания заявки
T_obsl_mean = T_obsl/length(ti); % Среднее время обслуживания заявки
T_being_mean= T_being/length(ti); % Среднее время нахождения заявки в системе

subplot(2,2,1);
semilogy(c, T_wait_mean, t, 'LineWidth', 2); xlabel('C'); ylabel('T_{WAIT} ( C )'); grid on; hold on
title('T_{WAIT} vs C');
subplot(2,2,2);
semilogy(U, T_wait_mean, t, 'LineWidth', 2); xlabel('U'); ylabel('T_{WAIT} ( U )'); grid on; hold on 
title('T_{WAIT} vs U');
subplot(2,2,3);
semilogy(c, L_wait_mean, t, 'LineWidth', 2); xlabel('C'); ylabel('L_{WAIT} ( C )'); grid on; hold on
title('L_{WAIT} vs C');
subplot(2,2,4);
semilogy(U, L_wait_mean, t, 'LineWidth', 2); xlabel('U'); ylabel('L_{WAIT} ( U )'); grid on; hold on
title('L_{WAIT} vs U');
figure(6);
% subplot(2,2,1);
% semilogy(c, T_obsl_mean, t, 'LineWidth', 2); xlabel('C'); ylabel('T_{SERV} ( C )'); grid on; hold on
% title('T_{SERV} vs C');
% subplot(2,2,2);
% semilogy(U, T_obsl_mean, t, 'LineWidth', 2); xlabel('U'); ylabel('T_{SERV} ( U )'); grid on; hold on
% title('T_{SERV} vs U');
subplot(2,1,1);
semilogy(c, T_being_mean, t, 'LineWidth', 2); xlabel('C'); ylabel('T_{BEING} ( C )'); grid on; hold on
title('T_{BEING} vs C');
subplot(2,1,2);
semilogy(U, T_being_mean, t, 'LineWidth', 2); xlabel('U'); ylabel('T_{BEING} ( U )'); grid on; hold on
title('T_{BEING} vs U');

end