%% 仿真1
% 简单LFM脉冲，脉压处理
close all;
clear;
clc;

%% 参数设置
f0 = 94e9;                  % 载频
T = 5e-6;                   % 脉宽
B = 20e6;                   % chirp信号频率调制带宽
C = 3e8;                    % 光速
K = B/T;                    % 调频率
Fs = 80e6;Ts = 1/Fs;        % 采样率和采样间隔
R = 100;                    % 目标距离
%% 产生回波
t = (1:round(T*Fs))/Fs;
% 目标回波
Srt = exp(1j*pi*K*(t-2*R/C).^2);
%% 用FFT和IFFT脉冲压缩处理
Nchirp = round(T/Ts);                   % 脉宽的点数
Srw = fft(Srt,Nchirp);                  % 雷达回波进行FFT
t0 = linspace(-T/2,T/2,Nchirp);
St = exp(1j*pi*K*t0.^2);                % 匹配滤波器时域响应
Sw = fft(St,Nchirp);                    % 匹配滤波器频域响应
Sot = fftshift(ifft(Srw.*conj(Sw)));    % 脉压后的信号
%========================================================================
% 绘图
figure;
subplot(3,1,1)
plot(t,real(Srt));axis tight;
xlabel('t/s');ylabel('幅度');title('雷达回波');
subplot(3,1,2)
plot(t0,real(St));axis tight;
xlabel('t/s');ylabel('幅度');title('匹配滤波器时域响应');
subplot(3,1,3)
plot(t*C/2,abs(Sot));
xlabel('距离/m');ylabel('幅度/dB');title('脉压后雷达回波');