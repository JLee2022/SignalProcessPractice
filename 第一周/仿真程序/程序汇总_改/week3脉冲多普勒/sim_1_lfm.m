%% 仿真1 LFM信号速度补偿
% 仿真目标运动对于脉冲压缩结果的影响，以及进行速度补偿
close all;
clear;
clc;

%% 参数设置
f0 = 94e9;                  % 载频
T = 5e-6;                   % 脉宽
B = 20e6;                   % chirp信号频率调制带宽
C = 3e8;                    % 光�??
K = B/T;                    % 调频�?
Fs = 80e6;Ts = 1/Fs;        % 采样率和采样间隔
R = 100;                    % 目标距离
V = 300;                     % 目标速度
lamda = C/f0;               % 波长
fd = 4*pi*V/lamda;          % 多普勒频�?
%% 产生回波
t = (1:round(T*Fs))/Fs;
% 目标回波
Srt = exp(1j*pi*K*(t-2*R/C).^2) .* exp(1j*fd*t);
%% 用FFT和IFFT脉冲压缩处理
Nchirp = round(T/Ts);                   % 脉宽的点�?
Srw = fft(Srt,Nchirp);                  % 雷达回波进行FFT
t0 = linspace(-T/2,T/2,Nchirp);

% 不进行�?�度补偿
St = exp(1j*pi*K*t0.^2);                % 匹配滤波器时域响�?
Sw = fft(St,Nchirp);                    % 匹配滤波器频域响�?
Sot = fftshift(ifft(Srw.*conj(Sw)));    % 脉压后的信号

% 进行速度补偿
St_v = exp(1j*pi*K*t0.^2) .* exp(1j*fd*t); % 匹配滤波器时域响�?
Sw_v = fft(St_v,Nchirp);                    % 匹配滤波器频域响�?
Sot_v = fftshift(ifft(Srw.*conj(Sw_v)));    % 脉压后的信号
%========================================================================
% 绘图
figure;
subplot(2,1,1)
plot(t*C/2,abs(Sot));
xlabel('距离/m');ylabel('幅度/dB');title('脉压后雷达回�?');
subplot(2,1,2)
plot(t*C/2,abs(Sot_v));
xlabel('距离/m');ylabel('幅度/dB');title('脉压后雷达回波（速度补偿�?');