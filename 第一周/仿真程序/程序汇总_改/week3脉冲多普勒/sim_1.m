%% 仿真1 LFM信号速度补偿
% 仿真目标运动对于去斜处理结果的影响，以及进行速度补偿
clc;
clear all;
close all;
%% 参数设置
f0 = 94e9;            % 中心频率
B = 20e6;             % 带宽
T = 5e-6;             % 脉宽
k = B/T;              % 调频率
Fs = 80e6;            % 采样率
R0 = 100;             % 目标距离
V0 = 100;
PRT = 50e-6;          % 脉冲重复时间
c = 3e8;

lamda = c/f0;         % 波长
fi = 4*pi*V0/lamda;   % 多普勒频率
N = round(PRT*Fs);    
fft_N = 2^nextpow2(N);
t = linspace(0,PRT,N);

%% 参考信号
% 目标无速度
Sref1=exp(2j*pi*f0*t).*exp(1j*pi*k*t.^2);
% 目标有速度
Sref2=exp(2j*pi*f0*t).*exp(1j*pi*k*t.^2).*exp(1j*fi*t);
%% 回波信号
% 目标有速度
St=exp(1j*pi*k*(t-2*R0/c).^2).*exp(2j*pi*f0*(t-2*R0/c)).*exp(1j*fi*t);
%% 混频信号
SSt1=St.*conj(Sref1);               %去斜后时域信号
spectrum1=fft(SSt1,fft_N);          %去斜后频域信号
f=Fs*(0:fft_N-1)/fft_N-Fs/2;        %从-Fs/2到Fs/2
f=f*c*T/2/B;                        %瞬时频率对应的距离
sf=exp(-1j*pi/k*f.^2);              %滤波器传输函数
SSt1=spectrum1.*sf;                 %从频域实现了压缩和去斜
SSt1=fftshift(SSt1);

%速度补偿后
SSt2=St.*conj(Sref2);               %去斜后时域信号
spectrum2=fft(SSt2,fft_N);          %去斜后频域信号
f=Fs*(0:fft_N-1)/fft_N-Fs/2;        %从-Fs/2到Fs/2
f=f*c*T/2/B;                        %瞬时频率对应的距离
sf=exp(-1j*pi/k*f.^2);              %滤波器传输函数
SSt2=spectrum2.*sf;                 %从频域实现了压缩和去斜
SSt2=fftshift(SSt2);

% 速度补偿对比结果作图
figure
subplot(2,1,1)
plot(abs(f),abs(SSt1));
xlabel('距离/m');title('速度补偿前');
grid on
subplot(2,1,2)
plot(abs(f),abs(SSt2));
xlabel('距离/m');title('速度补偿后');
grid on