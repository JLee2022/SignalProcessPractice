%% simu 1-1
close all
clear
%% ��������
N = 128;
PRT = 20e-6; % �����ظ����
fc = 20e6; % ����Ƶ��
fs = 4*fc; % �����ʣ�����������㻻Ϊ4������Ƶ�ʣ�
tau = 0.5*PRT; % ������Ϊһ���PRT

t = 0:1/fs:N*PRT; % һ֡�źŵĳ���
samp = length(t);
% ���������ź�
sig_rect = (square(2*pi*(1/PRT)*t, 25) + 1)/2;
sig_sin = sin(2*pi*fc*t);
sig_cos = cos(2*pi*fc*t);
sig_sin_rec = sig_rect.*sig_sin;
sig_cos_rec = sig_rect.*sig_cos;
sig_ref = cos(2*pi*fc*(t - 1.23e-6)); % �ο��ź�Ҫ��΢��һ���ӳ�

sig_I = sig_ref.*sig_cos;
sig_Q = sig_ref.*sig_sin; % IQͨ���ź�

% �����˲�������
Nf = 60; % �˲�������
f_pass = 10e6;
f_stop = 20e6;
wpass = 1;
wstop = 1;
dens = 20; % density factor
b = firpm(Nf, [0 f_pass f_stop fs/2]/(fs/2), [1 1 0 0], [wpass wstop], ...
            {dens});

% ��IQ�źŽ��е�ͨ�˲�����
fft_N = samp + Nf;
filter_I = ifft(fft(sig_I, fft_N).*fft(b, fft_N));
filter_Q = ifft(fft(sig_Q, fft_N).*fft(b, fft_N));
IQ = filter_I + 1j*filter_Q;

figure;
plot(linspace(-40, 40, samp), abs(fftshift(fft(sig_cos))));
xlabel('Ƶ��/MHz'); title('cos Signal');

figure;
subplot(2, 1, 1);
plot(linspace(-40, 40, 2*samp), abs(fftshift(fft(sig_I, 2*samp))));
xlabel('Ƶ��/MHz'); title('I Signal');
axis([-50,50,0,1000]);
subplot(2, 1, 2);
plot(linspace(-40, 40, 2*samp), abs(fftshift(fft(sig_Q, 2*samp))));
xlabel('Ƶ��/MHz'); title('Q Signal');
axis([-50,50,0,1000]);

figure;
subplot(2, 2, 1);
plot(linspace(-40, 40, fft_N), abs(fftshift(fft(filter_I, fft_N))));
xlabel('Ƶ��/MHz'); title('filter I Signal');
axis([-50,50,0,1000]);
subplot(2, 2, 2);
plot(linspace(-40, 40, fft_N), abs(fftshift(fft(filter_Q, fft_N))));
xlabel('Ƶ��/MHz'); title('Q Signal');
axis([-50,50,0,1000]);
subplot(2, 2, 3);
plot(linspace(-40, 40, fft_N), abs(fftshift(fft(b, fft_N))));
xlabel('Ƶ��/MHz'); title('filter');
subplot(2, 2, 4);
plot(linspace(-40, 40, fft_N), abs(fftshift(fft(IQ, fft_N))));
xlabel('Ƶ��/MHz'); title('IQ Signal');
axis([-50,50,0,1000]);


%% fir1����˲���: ����Ҫ����һ����ͨ�˲���
% N_fliter = 60; % �˲�������
% f_pass = [10e6, 20e6]; % ͨ�������
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


%% �˲������

