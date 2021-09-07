%% 仿真2 脉冲多普勒处理
clear;
close all;
clc;

%% 参数设置
C = 3e8;                %光速
Rinti = [100 500];      %初始距离
Vr = [50 100];          %速度
Np = 128;               %脉冲数量
prt = 20e-6;           %脉冲重复时间
f0 = 1e10;              %载频
lamda = C/f0;           %波长
Tw = 5e-6;              %脉冲时间
Bw = 20e6;              %带宽
K = Bw/Tw;              %调频率
Fs = 80e6;              %采样频率
Nsim = fix(Fs*prt);     %仿真点数
t_axis = (0:Nsim-1)*1/Fs;   %时间轴
r_axis = t_axis*C/2;        %距离轴


%% 脉冲多普勒
Sr = zeros(Np,Nsim);        %接收信号
Sp = zeros(Np,Nsim);        %脉压后信号
% 参考信号
Sref = rectpuls(t_axis-Tw/2,Tw).*exp(1i*pi*K.*(t_axis-Tw/2).^2);
Sref_fft = fft(Sref);
%% 步骤1:沿距离方向对每个脉冲进行脉冲压缩
% 生成回波并脉冲压缩
for n = 1:Np
    % 接收信号
    for i = 1:length(Rinti)
        tao = 2*(Rinti(i)-n*prt*Vr(i))/C;
        Sr1 = rectpuls(t_axis-tao-Tw/2,Tw).*exp(1i*pi*K.*(t_axis-tao-Tw/2).^2).*exp(-1i*2*pi*f0*tao);
        Sr1_fft = fft(Sr1);
        Sp1 = ifft(Sr1_fft.*conj(Sref_fft)); % 脉压后信号
        Sp(n,:) = Sp(n,:) + Sp1;
    end
end

% 画出脉压后信号的3D图
figure;
mesh(r_axis,1:Np,abs(Sp));
xlabel('距离[m]');ylabel('时间[prt]'); zlabel('幅度');
title('脉压后信号3D图');
grid on;

%% 步骤2: 沿慢时间做FFT
Spd = fft(Sp);
% 速度
fs_prt = 1/prt;                     %等效采样率
fd_axis = (0:Np-1)/Np*fs_prt;
v_axis = fd_axis*lamda/2;
v_resolution = 1/(Np*prt)*lamda/2;  %速度分辨率
% 画三维图
figure;
mesh(r_axis,v_axis,abs(Spd));
xlabel('距离[m]');ylabel('速度[m/s]');zlabel('幅度');
grid on;