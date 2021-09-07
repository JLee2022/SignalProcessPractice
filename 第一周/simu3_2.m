clear all;clc;close all;
fc=3e9;                 %�ز�Ƶ��
PRF=2000;       
Br=5e6;                 %����
fs=10*Br;               %����Ƶ��
Tp=5e-6;                %����
Kr=Br/Tp;               %Ƶ�ʱ仯��
c=3e8;                  %����
lamda=c/fc;             %����
Tr=1/PRF;               %�����ظ�����
N_mc=1.5/60*PRF;        %�������
t=0:1/fs:15*Tp+Tp;      %����ʱ��
N_r=length(t);          %��������
N_target=5;             %Ŀ�����
Rmax=c/2*15*Tp;                             %Ŀ�������루����Ӧ����1/2*c*Tr�����ǲ���ʱ�������˲�������ô��
R_t=Rmax*abs(rand(1,N_target))             %Ŀ��ľ��루��������Ŀ��ľ���һ����С��������ģ�
RCS_t=10*(exp(i*2*pi*rand(1,N_target)));    %Ŀ��RCS������Ϊ10����λ�ڣ�0,2pi��֮������ֲ�
Vmax=lamda*PRF/2;                           %Ŀ������ٶȣ������ٷ�Χ�����ڵ�һä��֮��
v=Vmax*((1+rand(1,N_target))/2)            %Ŀ���ٶȣ���������Ŀ����ٶ�һ����С�ڵ�һä�ٵģ���ÿһ��Ŀ�궼��һ���Լ����ٶȣ���Ӧһ������
%% ����Ŀ�����
sr=zeros(N_mc,N_r); %N_mc �������   N_r ��������
for i=1:N_mc
    ta=(i-1)*Tr;
    sri=0;%ÿһ�δ��ڲ�forѭ������֮��������Ϊ��һ������Ļز�������ŵ���һ������Ļز�
    %%�ڲ�forѭ����һ��Ŀ��һ��Ŀ�����о�����Ӧÿһ���ز���������ÿһ��Ŀ��ز�֮�����
    for k=1:N_target
        tao=2*(R_t(k)-v(k).*(ta+t))/c;
        srj=RCS_t(k).*rectpuls(t-tao-Tp/2,Tp).*exp(-1j*2*pi*fc*tao+1j*pi*Kr.*(t-tao-Tp/2).^2);
        sri=sri+srj;
    end
    %%���forѭ������ͬ�����壬��Ӧ��ta�ǲ�ֵͬ���ٴ���������ز�
    sr(i,:)=sri;
end

%% ����ѹ��ǰ�Ļز�
tm=(1:N_mc)/PRF;
R=c*t/2;
figure(1);
image(R,tm,255*abs(sr)/max(max(abs(sr))));
figure(2);
plot(t*c/2,abs(sr(1,:)));     %��ͼ����ֻ��ӳ�˵�һ������Ļز����

%% ����ѹ��
st=rectpuls(t-Tp/2,Tp).*exp(1i*pi*Kr*(t-Tp/2).^2);%�ο��ź� ʱ�� Ҳ����ƥ���˲�����ʱ��
stf=conj(fft(st));%ƥ���˲�����Ƶ������
for i=1:N_mc
    sr(i,:)=ifft(fft(sr(i,:)).*stf);  %�ֱ��ÿһ������ѹ�� Ƶ������ѹ��          
end
figure(3);
image(R,tm,255*abs(sr)/max(max(abs(sr)))) 
xlabel('����/m'); title('��ѹ��άĿ�������');
figure;
plot(t*c/2,abs(sr(1,:)))                       
             
sr=fft(sr,[],1);
V=linspace(0,PRF,50)*lamda/2;
figure;
image(R,V,255*abs(sr)/max(max(abs(sr))));
xlabel('����/m'); ylabel('�ٶȹ���ֵ'); title('������-������')
