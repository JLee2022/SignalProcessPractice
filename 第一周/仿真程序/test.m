%% 浠跨湡1 LFM淇″彿閫熷害琛ュ伩
% 浠跨湡鐩爣杩愬姩瀵逛簬鑴夊啿鍘嬬缉缁撴灉鐨勫奖鍝嶏紝浠ュ強杩涜閫熷害琛ュ伩
close all;
clear;
clc;

%% 鍙傛暟璁剧疆
f0 = 94e9;                  % 杞介
T = 5e-6;                   % 鑴夊
B = 20e6;                   % chirp淇″彿棰戠巼璋冨埗甯﹀
C = 3e8;                    % 鍏夐??
K = B/T;                    % 璋冮鐜?
Fs = 80e6;Ts = 1/Fs;        % 閲囨牱鐜囧拰閲囨牱闂撮殧
R = 100;                    % 鐩爣璺濈
V = 300;                     % 鐩爣閫熷害
lamda = C/f0;               % 娉㈤暱
fd = 4*pi*V/lamda;          % 澶氭櫘鍕掗鐜?
%% 浜х敓鍥炴尝
t = (1:round(T*Fs))/Fs;
% 鐩爣鍥炴尝
Srt = exp(1j*pi*K*(t-2*R/C).^2) .* exp(1j*fd*t);
%% 鐢‵FT鍜孖FFT鑴夊啿鍘嬬缉澶勭悊
Nchirp = round(T/Ts);                   % 鑴夊鐨勭偣鏁?
Srw = fft(Srt,Nchirp);                  % 闆疯揪鍥炴尝杩涜FFT
t0 = linspace(-T/2,T/2,Nchirp);

% 涓嶈繘琛岄?熷害琛ュ伩
St = exp(1j*pi*K*t0.^2);                % 鍖归厤婊ゆ尝鍣ㄦ椂鍩熷搷搴?
Sw = fft(St,Nchirp);                    % 鍖归厤婊ゆ尝鍣ㄩ鍩熷搷搴?
Sot = fftshift(ifft(Srw.*conj(Sw)));    % 鑴夊帇鍚庣殑淇″彿

% 杩涜閫熷害琛ュ伩
St_v = exp(1j*pi*K*t0.^2) .* exp(1j*fd*t); % 鍖归厤婊ゆ尝鍣ㄦ椂鍩熷搷搴?
Sw_v = fft(St_v,Nchirp);                    % 鍖归厤婊ゆ尝鍣ㄩ鍩熷搷搴?
Sot_v = fftshift(ifft(Srw.*conj(Sw_v)));    % 鑴夊帇鍚庣殑淇″彿
%========================================================================
% 缁樺浘
figure;
subplot(2,1,1)
plot(t*C/2,abs(Sot));
xlabel('璺濈/m');ylabel('骞呭害/dB');title('鑴夊帇鍚庨浄杈惧洖娉?');
subplot(2,1,2)
plot(t*C/2,abs(Sot_v));
xlabel('璺濈/m');ylabel('骞呭害/dB');title('鑴夊帇鍚庨浄杈惧洖娉紙閫熷害琛ュ伩锛?');










%% 仿真2
% close all;
% clear;
% clc;
% 
% %% 参数设置
% T = 5e-6;               % 脉宽
% B = 25e6;               % chirp信号频率调制带宽
% Rmin = 0;Rmax = 1000;   % 距离范围
% R = [100,110,300];      % 目标距离
% RCS = [1 1 3];          % 雷达散射面积
% C = 3e8;                % 光速
% K = B/T;                % 调频率
% Rwid = Rmax-Rmin;       % 距离接收窗
% Twid = 2*Rwid/C;        % 时间接收窗
% Fs = 100e6;Ts = 1/Fs;    % 采样率和采样间隔
% Nwid = ceil(Twid/Ts);   % 采样点接收窗，ceil朝正无穷方向取整
% %% 产生回波
% t = linspace(2*Rmin/C,2*Rmax/C,Nwid);               % 接收窗
% M = length(R);                                      % 目标数量
% td = ones(M,1)*t-2*R'/C*ones(1,Nwid); 
% Srt = RCS*(exp(1j*pi*K*td.^2).*(abs(td)<T/2));      % 目标回波
% %% 用FFT和IFFT脉冲压缩处理
% % 生成匹配滤波器
% Nchirp = round(T/Ts);                   % 脉宽的点数
% Nfft = 2^nextpow2(Nwid+Nwid-1);         % 最靠近括号内容的2的指数
% t0 = linspace(-T/2,T/2,Nchirp);
% St = exp(1j*pi*K*t0.^2);                % 匹配滤波器时域响应
% 
% % 匹配滤波处理
% St_fft = fft(St,Nfft);                              % 匹配滤波器频域响应
% St_fft_win = fft(St,Nfft) .* hamming(Nfft).';       % 加汉明窗的匹配滤波器频域响应
% Srt_fft = fft(Srt,Nfft);                            % 雷达回波进行FFT
% Sot = fftshift(ifft(Srt_fft.*conj(St_fft)));        % 脉压后的信号
% Sot_win = fftshift(ifft(Srt_fft.*conj(St_fft_win)));% 加窗脉压后的信号
% 
% % 根据接收窗选取信号
% N0=Nfft/2-Nchirp/2;
% Z=abs(Sot(N0:N0+Nwid-1));
% Z_win=abs(Sot_win(N0:N0+Nwid-1));
% 
% % 匹配滤波结果绘图
% figure;
% subplot(3,1,1)
% plot(t,real(Srt));axis tight;
% xlabel('t/s');ylabel('幅度');title('雷达回波');
% subplot(3,1,2)
% plot(t0,real(St))
% xlabel('t/s');ylabel('幅度');title('匹配滤波器时域响应');
% subplot(3,1,3)
% plot(t*C/2,Z);
% xlabel('距离/m');ylabel('幅度');title('脉压结果');
% 
% figure;
% subplot(2,1,1)
% plot(t*C/2,Z);
% xlabel('距离/m');ylabel('幅度');title('脉压结果');
% subplot(2,1,2)
% plot(t*C/2,Z_win);
% xlabel('距离/m');ylabel('幅度');title('脉压结果（加汉明窗）');













%% 数字下变频和滤波仿真2
% % 频率步进信号体制下，创建目标回波并对其进行数字下变频和数字滤波处理
% clear;clc;
% close all;
% 
% %% 参数设定
% c = 3e8;                    % 光速
% fs = 80e6;                  % 采样率
% tw = 0.1e-6;                % 脉宽
% f0 = 140e6;                 % 载频频率
% fc = 20e6 ;                 % 中心频率
% delt_f = 5e6;               % 步进频率
% Target_coordinate = [100];  % 目标距离
% RCS = [2];                  % 目标对应RCS
% PRT = 50e-6;                % 脉冲重复周期
% frame = 20e-5;              % 帧长度
% N = 128;                    % 一帧内脉冲个数
% 
% %% 建立回波信号
% % 采样点数为PRT*fs，表示在一个脉冲周期时间内有多少采样点
% x = zeros(N,round(PRT*fs)); 
% 
% % 生成一个脉冲长度范围内的时间坐标轴
% tt = (1:round(tw*fs))/fs;
% 
% for n = 1:N
%     % 建立每个PRT的回波信号，存储到回波矩阵中
%     for xi =1:length(Target_coordinate)  
%         y_tar = zeros(1,round(PRT*fs));
%         R = Target_coordinate(xi);
%         t_tar = round(2*R/c*fs)+1:round(2*R/c*fs)+round(tw*fs);
%         y_tar(t_tar) = RCS(xi) * exp(-1j*2*pi*(f0+(n-1)*delt_f)*(tt-2*R/c));
%         x(n,:) = x(n,:) + y_tar;
%     end 
% end
% 
% % 时域回波信号作图
% figure(1);
% plot((1:round(PRT*fs))/fs*c/2,abs(x(1,:))); 
% xlabel('距离/m');title('时域信号回波');
% 
% %% 下变频
% % 去除载频，目的是模拟信号经变频模块处理后，
% % 下变频之前所有PRT均保持相同的中心频率。
% y = zeros(N,round(PRT*fs));
% t = 1/fs:1/fs:PRT;
% fc = 20e6;     % 中心频率     
% for n = 1:N
%     y(n,:) = x(n,:) .* exp(-1j*2*pi*(-fc-f0-(n-1)*delt_f)*t);
% end
% 
% % 按照波门选取包含目标信息的部分进行处理
% data_range = 50:100;
% num_bomen = length(data_range);
% y_ddc_I = zeros(N,num_bomen);
% y_ddc_Q = zeros(N,num_bomen);
% for n = 1:N
%     y_ddc_I(n,:) = real(y(n,data_range)) .* cos(2*pi*fc*(0:num_bomen-1)/fs);
%     y_ddc_Q(n,:) = real(y(n,data_range)) .* sin(2*pi*fc*(0:num_bomen-1)/fs);
% end
% 
% % 下变频结果作图
% num = 2^13;
% figure;
% subplot(2,2,1);
% plot(linspace(-40,40,num), abs(fftshift(fft(x(1,:),num))));
% xlabel('频率/MHz');title('回波频域信号');
% subplot(2,2,2);
% plot(linspace(-40,40,num), abs(fftshift(fft(y(1,:),num))));
% xlabel('频率/MHz');title('去除载频后频域信号');
% subplot(2,2,3);
% plot(linspace(-40,40,num), abs(fftshift(fft(y_ddc_I(1,:),num))));
% xlabel('频率/MHz');title('I通道频域信号');
% subplot(2,2,4);
% plot(linspace(-40,40,num), abs(fftshift(fft(y_ddc_Q(1,:),num))));
% xlabel('频率/MHz');title('Q通道频域信号');
% 
% %% 滤波
% % 创建数字低通滤波器
% Nf=60;
% Fpass = 10e6;
% Fstop = 20e6; 
% Wpass = 1;    % Passband Weight
% Wstop = 1;    % Stopband Weight
% dens  = 20;   % Density Factor
% b  = firpm(Nf, [0 Fpass Fstop fs/2]/(fs/2), [1 1 0 0], [Wpass Wstop], ...
%            {dens});
% 
% % 进行滤波处理
% fft_N = num_bomen + Nf;
% data2_ddc_pulse_I = zeros(N,num_bomen);
% data2_ddc_pulse_Q = zeros(N,num_bomen);
% for n = 1:N
%     data2_filter_I = ifft(fft(y_ddc_I(n,:),fft_N).*fft(b,fft_N));
%     data2_filter_Q = ifft(fft(y_ddc_Q(n,:),fft_N).*fft(b,fft_N));
%     
%     data2_ddc_pulse_I(n,:) = data2_filter_I (1:num_bomen);
%     data2_ddc_pulse_Q(n,:) = data2_filter_Q (1:num_bomen);
% end
% data_ddc = data2_ddc_pulse_I + 1j*data2_ddc_pulse_Q;
% 
% % 数字低通滤波处理结果作图
% figure;
% subplot(2,2,1);
% plot(linspace(-40,40,fft_N), abs(fftshift(fft(y_ddc_I(1,:),fft_N))));
% xlabel('频率/MHz');title('I通道频域信号');
% subplot(2,2,2);
% plot(linspace(-40,40,fft_N), abs(fftshift(fft(b,fft_N))));
% xlabel('频率/MHz');title('低通滤波器频域信号');
% subplot(2,2,3);
% plot(linspace(-40,40,fft_N), abs(fftshift(fft(data_ddc(1,:),fft_N))));
% xlabel('频率/MHz');title('滤波结果频域信号');