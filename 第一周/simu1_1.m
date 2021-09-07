%% simu 1-1
close all
clear
%% 波形生成
N = 128;
PRT = 20e-6; % 脉冲重复间隔
fc = 20e6; % 中心频率
fs = 4*fc; % 采样率（方便仿真验算换为4倍中心频率）
tau = 0.5*PRT; % 脉冲宽度为一半的PRT

t = 0:1/fs:N*PRT; % 一帧信号的长度
samp = length(t);
% 生成脉冲信号
sig_rect = (square(2*pi*(1/PRT)*t, 25) + 1)/2;
sig_sin = sin(2*pi*fc*t);
sig_cos = cos(2*pi*fc*t);
sig_sin_rec = sig_rect.*sig_sin;
sig_cos_rec = sig_rect.*sig_cos;
sig_ref = cos(2*pi*fc*(t - 1.23e-6)); % 参考信号要稍微做一点延迟

sig_I = sig_ref.*sig_cos;
sig_Q = sig_ref.*sig_sin; % IQ通道信号

% 设置滤波器参数
Nf = 60; % 滤波器点数
f_pass = 10e6;
f_stop = 20e6;
wpass = 1;
wstop = 1;
dens = 20; % density factor
b = firpm(Nf, [0 f_pass f_stop fs/2]/(fs/2), [1 1 0 0], [wpass wstop], ...
            {dens});

% 对IQ信号进行低通滤波处理
fft_N = samp + Nf;
filter_I = ifft(fft(sig_I, fft_N).*fft(b, fft_N));
filter_Q = ifft(fft(sig_Q, fft_N).*fft(b, fft_N));
IQ = filter_I + 1j*filter_Q;

figure;
plot(linspace(-40, 40, samp), abs(fftshift(fft(sig_cos))));
xlabel('频率/MHz'); title('cos Signal');

figure;
subplot(2, 1, 1);
plot(linspace(-40, 40, 2*samp), abs(fftshift(fft(sig_I, 2*samp))));
xlabel('频率/MHz'); title('I Signal');
axis([-50,50,0,1000]);
subplot(2, 1, 2);
plot(linspace(-40, 40, 2*samp), abs(fftshift(fft(sig_Q, 2*samp))));
xlabel('频率/MHz'); title('Q Signal');
axis([-50,50,0,1000]);

figure;
subplot(2, 2, 1);
plot(linspace(-40, 40, fft_N), abs(fftshift(fft(filter_I, fft_N))));
xlabel('频率/MHz'); title('filter I Signal');
axis([-50,50,0,1000]);
subplot(2, 2, 2);
plot(linspace(-40, 40, fft_N), abs(fftshift(fft(filter_Q, fft_N))));
xlabel('频率/MHz'); title('Q Signal');
axis([-50,50,0,1000]);
subplot(2, 2, 3);
plot(linspace(-40, 40, fft_N), abs(fftshift(fft(b, fft_N))));
xlabel('频率/MHz'); title('filter');
subplot(2, 2, 4);
plot(linspace(-40, 40, fft_N), abs(fftshift(fft(IQ, fft_N))));
xlabel('频率/MHz'); title('IQ Signal');
axis([-50,50,0,1000]);


%% fir1设计滤波器: 仿真要求是一个带通滤波器
% N_fliter = 60; % 滤波器点数
% f_pass = [10e6, 20e6]; % 通带和阻带
% wn_pass = f_pass * 2 / fs;
% fli_bpass = fir1(N_fliter - 1, wn_pass, 'low');
% 
% re_bpass = 20*log(abs(fft(fli_bpass)))/log(10);
% xlabel_f = 0:(fs/length(re_bpass)):fs/2;
% figure
% yy = filter(fli_bpass, 1, sig);
% plot(yy)
% figure
% plot(abs(fft(yy, 256)))


%% 滤波器设计

