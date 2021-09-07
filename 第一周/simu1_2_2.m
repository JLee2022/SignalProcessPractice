%% simu1-2
clear all
close all
%% ���β���
c = 3e8; % ����
R = 100; % Ŀ�����
SNR = 10; % �����
tp = 10e-6; % ������
delta_f = 5e6; % ����Ƶ��
f_0 = 20e6; % ��ʼƵ��
fs = 400e6; % ������
N = 30; % ֡��������
PRT = 50e-6; PRI = 1/PRT; Tp = PRT; % �����ظ�����
t = 0:1/fs:(N)*PRT; % ����ʱ��(���۲��źŵ��ܳ���)
sig_t = zeros(length(t), 1)'; % �����ź�����
sig_r = zeros(length(t), 1)'; % �����ź�����
temp = zeros(length(t), 1)';
%% ���ɷ����ź�
for step = 1:N
    rec_freq = ((step - 1) * delta_f) + f_0; % ÿһ�θ��������崮��Ƶ��
    % ����ʱ��
    tau = (step - 1) * PRT;
    sig_t = sig_t + rectpuls(t - tau - tp/2, tp).*exp(j * 2 * pi * rec_freq * t);
end

%% ���ɻز��ź�
for step = 1:N
    rec_freq = ((step - 1) * delta_f) + f_0;
    % ����ʱ��
    tau = 2 * R / c + (step - 1) * PRT;
    sig_r = sig_r + rectpuls(t - tau - tp/2, tp).*exp(j * 2 * pi * rec_freq * (t-tau));
end

%% ��Ƶ
% ���ɲο��ź�
R_ref = 0; % ��ʹ�÷����źŽ��л�Ƶ
sig_ref = zeros(length(t), 1)';
for step = 1:N
    rec_freq = ((step - 1) * delta_f) + f_0;
    % ����ʱ��
    tau = 2 * R / c + (step - 1) * PRT;
    sig_ref= sig_ref + rectpuls(t - tau - tp/2, tp).*exp(j * 2 * pi * rec_freq * t);
end
figure
f = fs*((1:length(t)) - 1) / length(t);
plot(f, abs(fft(sig_r)))

% ��Ƶ����
sig_IQ = sig_r .* sig_ref;
figure
plot(f, abs(fft(sig_IQ)))
%% �˲�������
wp = 2e8; % ͨ����ֹƵ��
ws = 2.5e8; % �����ֹƵ��
ap = 0.1; as = 60;
wpp = wp/fs/2; wss = ws/fs/2;
[n, wn] = buttord(wpp, wss, ap, as);
[b, a] = butter(n, wn);
figure
freqz(b, a, 512, fs);
% �����˲����
y = filter(b, a, sig_IQ);
figure
plot(abs(fft(y)))
figure
plot(real(y))