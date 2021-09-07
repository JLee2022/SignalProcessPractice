%% 数字下变频和滤波仿真2
% 频率步进信号体制下，创建目标回波并对其进行数字下变频和数字滤波处理
clear;clc;
close all;

%% 参数设定
c = 3e8;                    % 光速
fs = 80e6;                  % 采样率
tw = 0.1e-6;                % 脉宽
f0 = 140e6;                 % 载频频率
fc = 20e6 ;                 % 中心频率
delt_f = 5e6;               % 步进频率
Target_coordinate = [100];  % 目标距离
RCS = [2];                  % 目标对应RCS
PRT = 50e-6;                % 脉冲重复周期
frame = 20e-5;              % 帧长度
N = 128;                    % 一帧内脉冲个数

%% 建立回波信号
% 采样点数为PRT*fs，表示在一个脉冲周期时间内有多少采样点
x = zeros(N,round(PRT*fs)); 

% 生成一个脉冲长度范围内的时间坐标轴
tt = (1:round(tw*fs))/fs;

for n = 1:N
    % 建立每个PRT的回波信号，存储到回波矩阵中
    for xi =1:length(Target_coordinate)  
        y_tar = zeros(1,round(PRT*fs));
        R = Target_coordinate(xi);
        t_tar = round(2*R/c*fs)+1:round(2*R/c*fs)+round(tw*fs);
        y_tar(t_tar) = RCS(xi) * exp(-1j*2*pi*(f0+(n-1)*delt_f)*(tt-2*R/c));
        x(n,:) = x(n,:) + y_tar;
    end 
end

% 时域回波信号作图
figure(1);
plot((1:round(PRT*fs))/fs*c/2,abs(x(1,:))); 
xlabel('距离/m');title('时域信号回波');

%% 下变频
% 去除载频，目的是模拟信号经变频模块处理后，
% 下变频之前所有PRT均保持相同的中心频率。
y = zeros(N,round(PRT*fs));
t = 1/fs:1/fs:PRT;
fc = 20e6;     % 中心频率     
for n = 1:N
    y(n,:) = x(n,:) .* exp(-1j*2*pi*(-fc-f0-(n-1)*delt_f)*t);
end

% 按照波门选取包含目标信息的部分进行处理
data_range = 50:100;
num_bomen = length(data_range);
y_ddc_I = zeros(N,num_bomen);
y_ddc_Q = zeros(N,num_bomen);
for n = 1:N
    y_ddc_I(n,:) = real(y(n,data_range)) .* cos(2*pi*fc*(0:num_bomen-1)/fs);
    y_ddc_Q(n,:) = real(y(n,data_range)) .* sin(2*pi*fc*(0:num_bomen-1)/fs);
end

% 下变频结果作图
num = 2^13;
figure;
subplot(2,2,1);
plot(linspace(-40,40,num), abs(fftshift(fft(x(1,:),num))));
xlabel('频率/MHz');title('回波频域信号');
subplot(2,2,2);
plot(linspace(-40,40,num), abs(fftshift(fft(y(1,:),num))));
xlabel('频率/MHz');title('去除载频后频域信号');
subplot(2,2,3);
plot(linspace(-40,40,num), abs(fftshift(fft(y_ddc_I(1,:),num))));
xlabel('频率/MHz');title('I通道频域信号');
subplot(2,2,4);
plot(linspace(-40,40,num), abs(fftshift(fft(y_ddc_Q(1,:),num))));
xlabel('频率/MHz');title('Q通道频域信号');

%% 滤波
% 创建数字低通滤波器
Nf=60;
Fpass = 10e6;
Fstop = 20e6; 
Wpass = 1;    % Passband Weight
Wstop = 1;    % Stopband Weight
dens  = 20;   % Density Factor
b  = firpm(Nf, [0 Fpass Fstop fs/2]/(fs/2), [1 1 0 0], [Wpass Wstop], ...
           {dens});

% 进行滤波处理
fft_N = num_bomen + Nf;
data2_ddc_pulse_I = zeros(N,num_bomen);
data2_ddc_pulse_Q = zeros(N,num_bomen);
for n = 1:N
    data2_filter_I = ifft(fft(y_ddc_I(n,:),fft_N).*fft(b,fft_N));
    data2_filter_Q = ifft(fft(y_ddc_Q(n,:),fft_N).*fft(b,fft_N));
    
    data2_ddc_pulse_I(n,:) = data2_filter_I (1:num_bomen);
    data2_ddc_pulse_Q(n,:) = data2_filter_Q (1:num_bomen);
end
data_ddc = data2_ddc_pulse_I + 1j*data2_ddc_pulse_Q;

% 数字低通滤波处理结果作图
figure;
subplot(2,2,1);
plot(linspace(-40,40,fft_N), abs(fftshift(fft(y_ddc_I(1,:),fft_N))));
xlabel('频率/MHz');title('I通道频域信号');
subplot(2,2,2);
plot(linspace(-40,40,fft_N), abs(fftshift(fft(b,fft_N))));
xlabel('频率/MHz');title('低通滤波器频域信号');
subplot(2,2,3);
plot(linspace(-40,40,fft_N), abs(fftshift(fft(data_ddc(1,:),fft_N))));
xlabel('频率/MHz');title('滤波结果频域信号');