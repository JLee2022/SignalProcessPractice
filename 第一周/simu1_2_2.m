%% simu1-2
clear all
close all
%% 波形参数
c = 3e8; % 光速
R = 100; % 目标距离
SNR = 10; % 信噪比
tp = 10e-6; % 脉冲宽度
delta_f = 5e6; % 步进频率
f_0 = 20e6; % 起始频率
fs = 400e6; % 采样率
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
    sig_t = sig_t + rectpuls(t - tau - tp/2, tp).*exp(j * 2 * pi * rec_freq * t);
end

%% 生成回波信号
for step = 1:N
    rec_freq = ((step - 1) * delta_f) + f_0;
    % 计算时延
    tau = 2 * R / c + (step - 1) * PRT;
    sig_r = sig_r + rectpuls(t - tau - tp/2, tp).*exp(j * 2 * pi * rec_freq * (t-tau));
end

%% 混频
% 生成参考信号
R_ref = 0; % 就使用发射信号进行混频
sig_ref = zeros(length(t), 1)';
for step = 1:N
    rec_freq = ((step - 1) * delta_f) + f_0;
    % 计算时延
    tau = 2 * R / c + (step - 1) * PRT;
    sig_ref= sig_ref + rectpuls(t - tau - tp/2, tp).*exp(j * 2 * pi * rec_freq * t);
end
figure
f = fs*((1:length(t)) - 1) / length(t);
plot(f, abs(fft(sig_r)))

% 混频处理
sig_IQ = sig_r .* sig_ref;
figure
plot(f, abs(fft(sig_IQ)))
%% 滤波器处理
wp = 2e8; % 通带截止频率
ws = 2.5e8; % 阻带截止频率
ap = 0.1; as = 60;
wpp = wp/fs/2; wss = ws/fs/2;
[n, wn] = buttord(wpp, wss, ap, as);
[b, a] = butter(n, wn);
figure
freqz(b, a, 512, fs);
% 计算滤波后的
y = filter(b, a, sig_IQ);
figure
plot(abs(fft(y)))
figure
plot(real(y))