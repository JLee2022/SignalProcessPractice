%pulstran����
%���ܣ��������з�����
%y=pulstran(t,d,'func')
%�ú�������һ����Ϊfunc��������������֮Ϊһ�����ڣ��Ӷ�����һ�������Ե���������(func�������Զ���)
%��pulstran�����ĺ����귶Χ������tָ����������d����ָ�������Ե�ƫ���������������ڵ����ĵ㣩������һ��func�����ᱻ����length(d)�Σ�
%�Ӷ�ʵ��һ�������������źŵĲ���
%y=pulstran(t,d,'func',p1,p2,...)
%���е�p1,p2,...Ϊ��Ҫת�͸�func�����Ķ����������ֵ�����˱���t֮�⣩

clear all
T=0:1/1E3:1;
D=0:1/4:1;
Y=pulstran(T,D,'rectpuls',0.1);%���������Եľ��������ź�
subplot(121)
plot(T,Y)
xlabel('t')
ylabel('h(t)')
grid on
axis([0 1 -0.1 1.1])
T=0:1/1E3:1;
D=0:1/3:1;
Y=pulstran(T,D,'tripuls',0.2,1);%���������Ե����ǲ��ź�
subplot(122)
plot(T,Y)
xlabel('t')
ylabel('h(t)')
grid on
axis([0 1 -0.1 1.1])
