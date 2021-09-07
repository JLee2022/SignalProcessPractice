%% simu1-2
%% 波形参数
close all
clear all
c = 3e8; % 光速
R = 1000; % 目标距离
SNR = 10; % 信噪比
tp = 10e-6; % 脉冲宽度
delta_f = 5e6; % 步进频率
f_0 = 20e6; % 起始频率
fs = 800e6; % 采样率
N = 30; % 帧内脉冲数
PRT = 50e-6; PRI = 1/PRT; Tp = PRT; % 脉冲重复周期
t = 0:1/fs:(N)*PRT; % 采样时间(即观测信号的总长度)
sig_t = zeros(length(t), 1)'; % 发射信号向量
sig_r = zeros(length(t), 1)'; % 接受信号向量
temp = zeros(length(t), 1)';
%% 生成发射信号
for step = 1:N
    rec_freq = ((step - 1) * delta_f) + f_0; % 每一次个发射脉冲串的频率
    % 计算时延
    tau = (step - 1) * PRT;
    sig_t = sig_t + rectpuls(t - tau - tp/2, tp).*cos(2 * pi * rec_freq * t);
end

%% 生成回波信号
for step = 1:N
    rec_freq = ((step - 1) * delta_f) + f_0;
    % 计算时延
    tau = 2 * R / c + (step - 1) * PRT;
    sig_r = sig_r + rectpuls(t - tau - tp/2, tp).*cos(2 * pi * rec_freq * (t - tau));
end
figure
plot(real(sig_r))
%% 混频
% 生成参考信号
R_ref = 0; % 就使用发射信号进行混频
sig_ref_I = zeros(length(t), 1)';
sig_ref_Q = zeros(length(t), 1)';
for step = 1:N
    rec_freq = ((step - 1) * delta_f) + f_0;
    % 计算时延
    tau = (step - 1) * PRT;
    sig_ref_I = sig_ref_I + rectpuls(t - tau - tp/2, tp).*cos(2 * pi * rec_freq * t);
end
for step = 1:N
    rec_freq = ((step - 1) * delta_f) + f_0;
    % 计算时延
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
% 混频处理
sig_I = sig_r*2.*sig_ref_I;
sig_Q = sig_r*2.*sig_ref_Q;
subplot(3, 1, 2);
plot(f, abs(fftshift(fft(sig_I))));
xlabel('freq/Hz'); title('I signal');
subplot(3, 1, 3);
plot(f, abs(fftshift(fft(sig_Q))));
xlabel('freq/Hz'); title('Q signal');
%% 滤波处理
f_stop = 20e6;
wn = (1/fs)*f_stop;
b1 = fir1(60, wn, 'low');
sig_I_lp = filter(b1, 1, sig_I);
sig_Q_lp = filter(b1, 1, sig_Q);
sig_IQ = sig_I_lp - 1j * sig_Q_lp;

figure
plot(f, abs(fftshift(fft(sig_IQ))));
xlabel('freq/Hz'); title('LPF IQ signal');


