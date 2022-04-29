function Potk = CMO_otkaz(vi, ti, c, U, t) % СМО с отказами
P = zeros(1, length(c));
parfor i = 1 : 1 : length(c)
    To = ti(1) + vi(1)/c(i); %Время обслуживания 1-го запроса
    for j = 2 : 1 : length(ti) % проверка доступа для текущей заявки
        if ti(j) < To % если время прихода j-го запроса меньше времени обслуживания текущего запроса
            P(i) = P(i) + 1; % то увеличиваем число отказов на один 
        else
            To = ti(j) + vi(j)/c(i); % принимаем запрос и записываем новое время, когда канал будет свободен 
        end
    end
end
Potk = P./length(ti); % матрица вероятностей отказа для разных С
subplot(2,1,1);
plot(c, Potk, t, 'LineWidth', 2); xlabel('C'); ylabel('P_{ОТКАЗ} ( C )'); grid on; title('P_{ОТКАЗ} vs C'); hold on
subplot(2,1,2);
plot(U, Potk, t, 'LineWidth', 2); xlabel('U'); ylabel('P_{ОТКАЗ} ( U )'); grid on; title('P_{ОТКАЗ} vs U'); hold on 
end
