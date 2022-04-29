clear all;
close all;
clc;

% ������� ��� ��. 6105. ������� 1. ���� 60. ������ + �������������
%% ���� �������� �������� ������
S = load('S.txt','%f'); %�������� � ������ S
S = S'; %��������.
S = S(1:1000);
V = load('V.txt','%f'); %�������� � ������ V
V = V'; %��������.
S_emp_sl = S(randperm(length(S))); %������ ������������

emp_dti_means = 1./S; % ������� �������� ���������� ����� ���������
emp_dti_means_sl = 1./S_emp_sl;
n = length(S); %����� �������

%% ��������� ������ �������� ������ ��� ������ ����������
poiss_emp_dti    = [ ];
poiss_emp_dti_sl = [ ];
poiss_emp_Vi     = [ ];
poiss_emp_Vi_sl  = [ ];

%% ������������ �������������
parfor j=1:n
   poiss_dti = exprnd(emp_dti_means(j),1,S(j)); %������������ ��������� ���������� ����� ���������
   poiss_dti_sl = exprnd(emp_dti_means_sl(j),1,S_emp_sl(j));
   k = randi(665370,1,S(j));
   k1 = randi(665370,1,S_emp_sl(j));
   poiss_Vi = V(k);% ��������� ������ ������� ��� ������ � 1 �������
   poiss_Vi_sl = V(k1);
   poiss_emp_dti = [poiss_emp_dti poiss_dti];% ������� dti
   poiss_emp_dti_sl = [poiss_emp_dti_sl  poiss_dti_sl];
   poiss_emp_Vi = [poiss_emp_Vi poiss_Vi]; % ������� emp_Vi
   poiss_emp_Vi_sl = [poiss_emp_Vi_sl poiss_Vi_sl];
end
emp_ti = cumsum(poiss_emp_dti);
emp_ti_sl = cumsum(poiss_emp_dti_sl);
emp_ti = emp_ti-min(emp_ti); % �������� ��������
emp_ti_sl = emp_ti_sl - min(emp_ti_sl); % ��� �������� ������������ ������������������ �� ��������
emp_Vi = poiss_emp_Vi; % �������� ������
emp_Vi_sl = poiss_emp_Vi_sl;

figure % ������ ����������� ����������� ��������
stem(emp_ti(1:5000),log10(emp_Vi(1:5000)),'k','LineWidth',1);
grid minor;
title('����������� ����������� ����������� ��������');
xlabel('ti');
ylabel('log_{10}( Vi(ti) )');

c0 = sum(emp_Vi)/sum(poiss_emp_dti); % ������� ���������� ����������� 
c = 1.1*c0 : c0/2 : 10*c0; % ������� ���������� ����������� ������
U = c0./c; % ������� ����������� ������������� ������� ������


%% ��� ���� M/M/1
ne = length(poiss_emp_dti);
exp_dti = mean(poiss_emp_dti);
exp_dti = exprnd(exp_dti,1,ne);
exp_ti = cumsum(exp_dti);
exp_Vi = exprnd(mean(emp_Vi),1,ne);

%% ��� � �������� �������������� ����. ���������� (������)
PH = gpfit(poiss_emp_dti); % ����������� ���������� ������ (��������������.)
theta = PH(2)/PH(1);
GP_dtij = gprnd(PH(1),PH(2),theta,1,ne);
GP_ti = cumsum(GP_dtij);
GP_ti = (GP_ti - mean(GP_ti)) / std(GP_ti) * std(emp_ti) + mean(emp_ti); %����������

PH = gpfit(S); % ����������� ���������� �������� S ��� ������
GP_int = gprnd(PH(1),PH(2),PH(2)/PH(1),1,n);
GP_int = (GP_int - mean(GP_int)) / std(GP_int) * std(S) + mean(S); %����������

%% ��� � �������� �������������� ������� (�������������)
LoGNormn_prmtr = lognfit(emp_Vi); %����������� ��������� �������������
LogN_Vi = lognrnd(LoGNormn_prmtr(1),LoGNormn_prmtr(2),1,ne);
LogN_Vi = (LogN_Vi - mean(LogN_Vi)) / std(LogN_Vi) * std(emp_Vi) + mean(emp_Vi);

%% �������� ��������-������
% ������ ���������� ������
N = input('������� ���������� ������� ����������������� �������� = ');
H = Herst1(S-mean(S),N);
%disp(H);
H = H(2,1);

% ���������� ��������� ������� DFA
S_mod2 = alg_SHSH(GP_int, H, n, N);

%% ������������ ������ � ����������� � �������������� ������
GP_dti = [ ];
for i=1:length(S_mod2)
    GP_dti1 = gprnd(1/S_mod2(i), PH(2), theta, 1, S_emp_sl(i));
    GP_dti  = [GP_dti GP_dti1];
end
GP_tiS = cumsum(GP_dti);
GP_tiS = GP_tiS - min(GP_tiS);
GP_tiS = (GP_tiS - mean(GP_tiS)) / std(GP_tiS) * std(emp_ti) + mean(emp_ti); % ����������

%% ������ ��� ��� � ����������
ro_p = std(poiss_emp_dti)/mean(poiss_emp_dti); % ����������� �������� dti
ro_v = std(emp_Vi)/mean(emp_Vi);% ����������� �������� Vi
T_wait_mean_teor = (mean(emp_Vi)/c0)*((U.^2)./(1-U)).*((ro_p^2)+(ro_v^2))./2; % ������� ����� ��������
L_wait_mean_teor = (1/mean(poiss_emp_dti))*T_wait_mean_teor; % ������� ����� �������

%% ��� � ����������
figure(5);
figure(6)
t ='r--'; [T L] = CMO_ozhid(emp_Vi, emp_ti, c, U, t);
t ='g--'; [T L] = CMO_ozhid(exp_Vi, exp_ti, c, U, t);
t ='b--'; [T L] = CMO_ozhid(emp_Vi, GP_ti, c, U, t);
t ='k--'; [T L] = CMO_ozhid(emp_Vi, GP_tiS, c, U, t);
t ='m--'; [T L] = CMO_ozhid(LogN_Vi, emp_ti, c, U, t);
t ='c--'; [T L] = CMO_ozhid(emp_Vi_sl, emp_ti_sl, c, U, t);
%t ='y--o'; [T L] = CMO_ozhid(puass_Vi, emp_ti, c, U, t);
hold on;
legend('emp/emp/1', 'exp/exp/1', 'Pareto/emp/1', 'Pareto ��/emp/1', 'emp/LogNormal/1', 'theor.');

figure(5);
subplot(2,2,1);
semilogy(c, T_wait_mean_teor, 'y--o', 'LineWidth', 2); xlabel('C'); ylabel('T_{WAIT} ( C )'); grid minor; hold on
subplot(2,2,2);
semilogy(U, T_wait_mean_teor, 'y--o', 'LineWidth', 2); xlabel('U'); ylabel('T_{WAIT} ( U )'); grid minor; hold on 
subplot(2,2,3);
semilogy(c, L_wait_mean_teor, 'y--o', 'LineWidth', 2); xlabel('C'); ylabel('L_{WAIT} ( C )'); grid minor; hold on
subplot(2,2,4);
semilogy(U, L_wait_mean_teor, 'y--o', 'LineWidth', 2); xlabel('U'); ylabel('L_{WAIT} ( U )'); grid minor; hold on 

legend('emp/emp/1', 'exp/exp/1', 'Pareto/emp/1', 'Pareto ��/emp/1', 'emp/LogNorm/1', 'sl_{emp}/sl_{emp}/1', 'theor.');
figure(6);
legend('emp/emp/1', 'exp/exp/1', 'Pareto/emp/1', 'Pareto ��/emp/1', 'emp/LogNorm/1', 'sl_{emp}/sl_{emp}/1');

%% ��� � �������� 
figure('name', '�������������� �������������� ��� ��� � ��������')
t = 'r--';  emp_Potk = CMO_otkaz(emp_Vi, emp_ti, c,U,t); 
t = 'g--'; exp_Potk = CMO_otkaz(exp_Vi, exp_ti, c,U,t); 
t = 'b--'; LogN_Potk = CMO_otkaz(emp_Vi, GP_ti, c,U,t);    
t = 'k--'; LogN_SHSH_Potk = CMO_otkaz(emp_Vi, GP_tiS, c,U,t);
t = 'm--'; Bino_Potk = CMO_otkaz(LogN_Vi, emp_ti, c,U,t);
t = 'c--'; sl_Potk = CMO_otkaz(emp_Vi_sl, emp_ti_sl, c,U,t);
%t = 'y--o'; PUASS_Potk = CMO_otkaz(puass_Vi, emp_ti, c,U,t);
legend('emp/emp/1', 'exp/exp/1', 'Pareto/emp/1', 'Pareto ��/emp/1', 'emp/LogNorm/1', 'sl_{emp}/sl_{emp}/1');

% �������������� ��������������
%% ��������� ���������
q='r-'; 
g=Statistic(poiss_emp_dti,    1, q, '����. �������������� ��� emp dti',0);
g=Statistic(GP_dtij,          1, q, '����. �������������� ��� Pareto dti',0);
g=Statistic(GP_dti,       'gp' , q, '����. �������������� ��� Pareto dti (��)',0);
g=Statistic(exp_dti,      'exp', q, '����. �������������� ��� exp dti',0);
g=Statistic(poiss_emp_dti_sl, 1, q, '����. �������������� ��� rnd dti',0);

figure('name','������� ������������� ����. ����������'); hold on; grid minor;
q='r-';  g=build_cdf(poiss_emp_dti,q);
q='b--'; g=build_cdf(GP_dtij,q);
q='g--'; g=build_cdf(GP_dti,q);
q='c--'; g=build_cdf(exp_dti,q);
q='k--'; g=build_cdf(poiss_emp_dti_sl,q);
legend('�������','�������-�����','������','EXP');
xlabel('{\it X}'); ylabel('{\it F ( X )}'); title('������� ������������� ����. ����������');

%% ������ ������
q='r-'; 
g=Statistic(emp_Vi,     1, q,    '����. �������������� ��� emp Vi',1);
g=Statistic(LogN_Vi,    1, q,    '����. �������������� ��� LogNorm Vi',1);
g=Statistic(emp_Vi_sl,  1, q,    '����. �������������� ��� rnd Vi',1);
g=Statistic(exp_Vi, 'exp', q,    '����. �������������� ��� exp Vi',0);


figure('name','������� ������������� ������� ������'); hold on; grid minor;
q='r-';  g=build_cdf(emp_Vi,q);
q='b--'; g=build_cdf(LogN_Vi,q);
q='g--'; g=build_cdf(emp_Vi_sl,q);
q='k--'; g=build_cdf(exp_Vi,q);
legend('��������','�������������','�������-�����','M/M/1');
xlabel('{\it X}'); ylabel('{\it F ( X )}'); title('������� ������������� ������ ������');
xlim([0 5000]);

%% �������������
q='r-';
g=Statistic(S,        1, q, '����. �������������� ��� �������. ������.',0);
g=Statistic(S_emp_sl, 1, q, '����. �������������� ��� ����. ������.',0);
g=Statistic(GP_int,   1, q, '����. �������������� ��� Pareto ������.',0);
g=Statistic(S_mod2,   1, q, '����. �������������� ��� Pareto ������. � ��',0);


figure('name','������� ������������� ������.'); hold on; grid minor;
q='r-';  g=build_cdf(S,q);
q='b--'; g=build_cdf(S_emp_sl,q);
q='g--'; g=build_cdf(GP_int,q);
%q='k--'; g=build_cdf(S_mod2,q);
legend('��������','�������-�����','������');
xlabel('{\it X}'); ylabel('{\it F ( X )}'); title('������� ������������� ���������.');

disp(['� ���-�� ��������, � ��� - �� �����'])