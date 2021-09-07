%% simu1-2
%% ���β���
close all
clear all
c = 3e8; % ����
R = 1000; % Ŀ�����
SNR = 10; % �����
tp = 10e-6; % ������
delta_f = 5e6; % ����Ƶ��
f_0 = 20e6; % ��ʼƵ��
fs = 800e6; % ������
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
    sig_t = sig_t + rectpuls(t - tau - tp/2, tp).*cos(2 * pi * rec_freq * t);
end

%% ���ɻز��ź�
for step = 1:N
    rec_freq = ((step - 1) * delta_f) + f_0;
    % ����ʱ��
    tau = 2 * R / c + (step - 1) * PRT;
    sig_r = sig_r + rectpuls(t - tau - tp/2, tp).*cos(2 * pi * rec_freq * (t - tau));
end
figure
plot(real(sig_r))
%% ��Ƶ
% ���ɲο��ź�
R_ref = 0; % ��ʹ�÷����źŽ��л�Ƶ
sig_ref_I = zeros(length(t), 1)';
sig_ref_Q = zeros(length(t), 1)';
for step = 1:N
    rec_freq = ((step - 1) * delta_f) + f_0;
    % ����ʱ��
    tau = (step - 1) * PRT;
    sig_ref_I = sig_ref_I + rectpuls(t - tau - tp/2, tp).*cos(2 * pi * rec_freq * t);
end
for step = 1:N
    rec_freq = ((step - 1) * delta_f) + f_0;
    % ����ʱ��
    tau = (step - 1) * PRT;
    sig_ref_Q = sig_ref_Q + rectpuls(t - tau - tp/2, tp).*sin(2 * pi * rec_freq * t);
end
f = fs*((-length(t)/2:length(t)/2-1)) / length(t);
% figure
% plot(sig_r)
figure
subplot(3, 1, 1)
plot(f, abs(fftshift(fft(sig_r))));
xlabel('freq/Hz'); title('echo signal');
% ��Ƶ����
sig_I = sig_r*2.*sig_ref_I;
sig_Q = sig_r*2.*sig_ref_Q;
subplot(3, 1, 2);
plot(f, abs(fftshift(fft(sig_I))));
xlabel('freq/Hz'); title('I signal');
subplot(3, 1, 3);
plot(f, abs(fftshift(fft(sig_Q))));
xlabel('freq/Hz'); title('Q signal');
%% �˲�����
f_stop = 20e6;
wn = (1/fs)*f_stop;
b1 = fir1(60, wn, 'low');
sig_I_lp = filter(b1, 1, sig_I);
sig_Q_lp = filter(b1, 1, sig_Q);
sig_IQ = sig_I_lp - 1j * sig_Q_lp;

figure
plot(f, abs(fftshift(fft(sig_IQ))));
xlabel('freq/Hz'); title('LPF IQ signal');


