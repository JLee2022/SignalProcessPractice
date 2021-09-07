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










%% ����2
% close all;
% clear;
% clc;
% 
% %% ��������
% T = 5e-6;               % ����
% B = 25e6;               % chirp�ź�Ƶ�ʵ��ƴ���
% Rmin = 0;Rmax = 1000;   % ���뷶Χ
% R = [100,110,300];      % Ŀ�����
% RCS = [1 1 3];          % �״�ɢ�����
% C = 3e8;                % ����
% K = B/T;                % ��Ƶ��
% Rwid = Rmax-Rmin;       % ������մ�
% Twid = 2*Rwid/C;        % ʱ����մ�
% Fs = 100e6;Ts = 1/Fs;    % �����ʺͲ������
% Nwid = ceil(Twid/Ts);   % ��������մ���ceil���������ȡ��
% %% �����ز�
% t = linspace(2*Rmin/C,2*Rmax/C,Nwid);               % ���մ�
% M = length(R);                                      % Ŀ������
% td = ones(M,1)*t-2*R'/C*ones(1,Nwid); 
% Srt = RCS*(exp(1j*pi*K*td.^2).*(abs(td)<T/2));      % Ŀ��ز�
% %% ��FFT��IFFT����ѹ������
% % ����ƥ���˲���
% Nchirp = round(T/Ts);                   % ����ĵ���
% Nfft = 2^nextpow2(Nwid+Nwid-1);         % ����������ݵ�2��ָ��
% t0 = linspace(-T/2,T/2,Nchirp);
% St = exp(1j*pi*K*t0.^2);                % ƥ���˲���ʱ����Ӧ
% 
% % ƥ���˲�����
% St_fft = fft(St,Nfft);                              % ƥ���˲���Ƶ����Ӧ
% St_fft_win = fft(St,Nfft) .* hamming(Nfft).';       % �Ӻ�������ƥ���˲���Ƶ����Ӧ
% Srt_fft = fft(Srt,Nfft);                            % �״�ز�����FFT
% Sot = fftshift(ifft(Srt_fft.*conj(St_fft)));        % ��ѹ����ź�
% Sot_win = fftshift(ifft(Srt_fft.*conj(St_fft_win)));% �Ӵ���ѹ����ź�
% 
% % ���ݽ��մ�ѡȡ�ź�
% N0=Nfft/2-Nchirp/2;
% Z=abs(Sot(N0:N0+Nwid-1));
% Z_win=abs(Sot_win(N0:N0+Nwid-1));
% 
% % ƥ���˲������ͼ
% figure;
% subplot(3,1,1)
% plot(t,real(Srt));axis tight;
% xlabel('t/s');ylabel('����');title('�״�ز�');
% subplot(3,1,2)
% plot(t0,real(St))
% xlabel('t/s');ylabel('����');title('ƥ���˲���ʱ����Ӧ');
% subplot(3,1,3)
% plot(t*C/2,Z);
% xlabel('����/m');ylabel('����');title('��ѹ���');
% 
% figure;
% subplot(2,1,1)
% plot(t*C/2,Z);
% xlabel('����/m');ylabel('����');title('��ѹ���');
% subplot(2,1,2)
% plot(t*C/2,Z_win);
% xlabel('����/m');ylabel('����');title('��ѹ������Ӻ�������');













%% �����±�Ƶ���˲�����2
% % Ƶ�ʲ����ź������£�����Ŀ��ز���������������±�Ƶ�������˲�����
% clear;clc;
% close all;
% 
% %% �����趨
% c = 3e8;                    % ����
% fs = 80e6;                  % ������
% tw = 0.1e-6;                % ����
% f0 = 140e6;                 % ��ƵƵ��
% fc = 20e6 ;                 % ����Ƶ��
% delt_f = 5e6;               % ����Ƶ��
% Target_coordinate = [100];  % Ŀ�����
% RCS = [2];                  % Ŀ���ӦRCS
% PRT = 50e-6;                % �����ظ�����
% frame = 20e-5;              % ֡����
% N = 128;                    % һ֡���������
% 
% %% �����ز��ź�
% % ��������ΪPRT*fs����ʾ��һ����������ʱ�����ж��ٲ�����
% x = zeros(N,round(PRT*fs)); 
% 
% % ����һ�����峤�ȷ�Χ�ڵ�ʱ��������
% tt = (1:round(tw*fs))/fs;
% 
% for n = 1:N
%     % ����ÿ��PRT�Ļز��źţ��洢���ز�������
%     for xi =1:length(Target_coordinate)  
%         y_tar = zeros(1,round(PRT*fs));
%         R = Target_coordinate(xi);
%         t_tar = round(2*R/c*fs)+1:round(2*R/c*fs)+round(tw*fs);
%         y_tar(t_tar) = RCS(xi) * exp(-1j*2*pi*(f0+(n-1)*delt_f)*(tt-2*R/c));
%         x(n,:) = x(n,:) + y_tar;
%     end 
% end
% 
% % ʱ��ز��ź���ͼ
% figure(1);
% plot((1:round(PRT*fs))/fs*c/2,abs(x(1,:))); 
% xlabel('����/m');title('ʱ���źŻز�');
% 
% %% �±�Ƶ
% % ȥ����Ƶ��Ŀ����ģ���źž���Ƶģ�鴦���
% % �±�Ƶ֮ǰ����PRT��������ͬ������Ƶ�ʡ�
% y = zeros(N,round(PRT*fs));
% t = 1/fs:1/fs:PRT;
% fc = 20e6;     % ����Ƶ��     
% for n = 1:N
%     y(n,:) = x(n,:) .* exp(-1j*2*pi*(-fc-f0-(n-1)*delt_f)*t);
% end
% 
% % ���ղ���ѡȡ����Ŀ����Ϣ�Ĳ��ֽ��д���
% data_range = 50:100;
% num_bomen = length(data_range);
% y_ddc_I = zeros(N,num_bomen);
% y_ddc_Q = zeros(N,num_bomen);
% for n = 1:N
%     y_ddc_I(n,:) = real(y(n,data_range)) .* cos(2*pi*fc*(0:num_bomen-1)/fs);
%     y_ddc_Q(n,:) = real(y(n,data_range)) .* sin(2*pi*fc*(0:num_bomen-1)/fs);
% end
% 
% % �±�Ƶ�����ͼ
% num = 2^13;
% figure;
% subplot(2,2,1);
% plot(linspace(-40,40,num), abs(fftshift(fft(x(1,:),num))));
% xlabel('Ƶ��/MHz');title('�ز�Ƶ���ź�');
% subplot(2,2,2);
% plot(linspace(-40,40,num), abs(fftshift(fft(y(1,:),num))));
% xlabel('Ƶ��/MHz');title('ȥ����Ƶ��Ƶ���ź�');
% subplot(2,2,3);
% plot(linspace(-40,40,num), abs(fftshift(fft(y_ddc_I(1,:),num))));
% xlabel('Ƶ��/MHz');title('Iͨ��Ƶ���ź�');
% subplot(2,2,4);
% plot(linspace(-40,40,num), abs(fftshift(fft(y_ddc_Q(1,:),num))));
% xlabel('Ƶ��/MHz');title('Qͨ��Ƶ���ź�');
% 
% %% �˲�
% % �������ֵ�ͨ�˲���
% Nf=60;
% Fpass = 10e6;
% Fstop = 20e6; 
% Wpass = 1;    % Passband Weight
% Wstop = 1;    % Stopband Weight
% dens  = 20;   % Density Factor
% b  = firpm(Nf, [0 Fpass Fstop fs/2]/(fs/2), [1 1 0 0], [Wpass Wstop], ...
%            {dens});
% 
% % �����˲�����
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
% % ���ֵ�ͨ�˲���������ͼ
% figure;
% subplot(2,2,1);
% plot(linspace(-40,40,fft_N), abs(fftshift(fft(y_ddc_I(1,:),fft_N))));
% xlabel('Ƶ��/MHz');title('Iͨ��Ƶ���ź�');
% subplot(2,2,2);
% plot(linspace(-40,40,fft_N), abs(fftshift(fft(b,fft_N))));
% xlabel('Ƶ��/MHz');title('��ͨ�˲���Ƶ���ź�');
% subplot(2,2,3);
% plot(linspace(-40,40,fft_N), abs(fftshift(fft(data_ddc(1,:),fft_N))));
% xlabel('Ƶ��/MHz');title('�˲����Ƶ���ź�');