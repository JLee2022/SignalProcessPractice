%% 仿真3 去斜处理
% 仿真大宽带LFM波去斜处理流程
clc;
clear all;
close all;
%% 参数设置
f0 = 94e9;            %中心频率
B = 20e6;             %带宽
T = 5e-6;            %脉宽
k = B/T;              %调频率
Fs = 80e6;            %采样率
R = [100,110,150];%目标距离
PRT = 30e-6;          %脉冲重复时间
c = 3e8;

N = round(PRT*Fs);
fft_N = 2^nextpow2(N);
t = linspace(0,PRT,N);

%% 参考信号
Sref=exp(2j*pi*f0*t).*exp(1j*pi*k*t.^2);
%% 回波信号
St = zeros(1,N);
for i = 1:length(R)
    St_i = exp(1j*pi*k*(t-2*R(i)/c).^2).*exp(2j*pi*f0*(t-2*R(i)/c));
    St = St + St_i;
end
%% 混频信号
SSt = St.*conj(Sref);             % 去除载频
f = Fs*(0:fft_N-1)/fft_N;    
f = f*c*T/2/B;                    % 瞬时频率转换为对应的距离
% sf = exp(-1j*pi*k*f.^2);          % 滤波器传输函数
SSt = fliplr(fft(SSt,fft_N));         % 从频域实现了压缩和去斜

% 去斜结果作图
figure;
subplot(2,1,1);
plot(f,db(abs(SSt)/max(SSt)))
xlabel('距离/m');title('去斜处理前');
grid on
subplot(2,1,2);
plot(f,abs(SSt))
xlabel('距离/m');title('去斜处理后');
grid on