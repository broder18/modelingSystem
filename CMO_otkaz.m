function Potk = CMO_otkaz(vi, ti, c, U, t) % ��� � ��������
P = zeros(1, length(c));
parfor i = 1 : 1 : length(c)
    To = ti(1) + vi(1)/c(i); %����� ������������ 1-�� �������
    for j = 2 : 1 : length(ti) % �������� ������� ��� ������� ������
        if ti(j) < To % ���� ����� ������� j-�� ������� ������ ������� ������������ �������� �������
            P(i) = P(i) + 1; % �� ����������� ����� ������� �� ���� 
        else
            To = ti(j) + vi(j)/c(i); % ��������� ������ � ���������� ����� �����, ����� ����� ����� �������� 
        end
    end
end
Potk = P./length(ti); % ������� ������������ ������ ��� ������ �
subplot(2,1,1);
plot(c, Potk, t, 'LineWidth', 2); xlabel('C'); ylabel('P_{�����} ( C )'); grid on; title('P_{�����} vs C'); hold on
subplot(2,1,2);
plot(U, Potk, t, 'LineWidth', 2); xlabel('U'); ylabel('P_{�����} ( U )'); grid on; title('P_{�����} vs U'); hold on 
end
