%% 仿真2
close all;
clear;
clc;

%% 参数设置
T = 5e-6;               % 脉宽
B = 25e6;               % chirp信号频率调制带宽
Rmin = 0;Rmax = 1000;   % 距离范围
R = [100,110,150];      % 目标距离
RCS = [1 1 1];          % 雷达散射面积
C = 3e8;                % 光速
K = B/T;                % 调频率
Rwid = Rmax-Rmin;       % 距离接收窗
Twid = 2*Rwid/C;        % 时间接收窗
Fs = 80e6;Ts = 1/Fs;    % 采样率和采样间隔
Nwid = ceil(Twid/Ts);   % 采样点接收窗，ceil朝正无穷方向取整
%% 产生回波
t = linspace(2*Rmin/C,2*Rmax/C,Nwid);               % 接收窗
M = length(R);                                      % 目标数量
td = ones(M,1)*t-2*R'/C*ones(1,Nwid); 
Srt = RCS*(exp(1j*pi*K*td.^2).*(abs(td)<T/2));      % 目标回波
%% 用FFT和IFFT脉冲压缩处理
% 生成匹配滤波器
Nchirp = round(T/Ts);                   % 脉宽的点数
Nfft = 2^nextpow2(Nwid+Nwid-1);         % 最靠近括号内容的2的指数
t0 = linspace(-T/2,T/2,Nchirp);
St = exp(1j*pi*K*t0.^2);                % 匹配滤波器时域响应

% 匹配滤波处理
St_fft = fft(St,Nfft);                              % 匹配滤波器频域响应
St_fft_win = fft(St,Nfft) .* hamming(Nfft).';       % 加汉明窗的匹配滤波器频域响应
Srt_fft = fft(Srt,Nfft);                            % 雷达回波进行FFT
Sot = fftshift(ifft(Srt_fft.*conj(St_fft)));        % 脉压后的信号
Sot_win = fftshift(ifft(Srt_fft.*conj(St_fft_win)));% 加窗脉压后的信号

% 根据接收窗选取信号
N0=Nfft/2-Nchirp/2;
Z=abs(Sot(N0:N0+Nwid-1));
Z_win=abs(Sot_win(N0:N0+Nwid-1));

% 匹配滤波结果绘图
figure;
subplot(3,1,1)
plot(t,real(Srt));axis tight;
xlabel('t/s');ylabel('幅度');title('雷达回波');
subplot(3,1,2)
plot(t0,real(St))
xlabel('t/s');ylabel('幅度');title('匹配滤波器时域响应');
subplot(3,1,3)
plot(t*C/2,Z);
xlabel('距离/m');ylabel('幅度');title('脉压结果');

figure;
subplot(2,1,1)
plot(t*C/2,Z);
xlabel('距离/m');ylabel('幅度');title('脉压结果');
subplot(2,1,2)
plot(t*C/2,Z_win);
xlabel('距离/m');ylabel('幅度');title('脉压结果（加汉明窗）');