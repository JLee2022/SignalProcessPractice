%pulstran函数
%功能：脉冲序列发生器
%y=pulstran(t,d,'func')
%该函数基于一个名为func的连续函数并以之为一个周期，从而产生一串周期性的连续函数(func函数可自定义)
%该pulstran函数的横坐标范围由向量t指定，而向量d用于指定周期性的偏移量（即各个周期的中心点），这样一个func函数会被计算length(d)次，
%从而实现一个周期性脉冲信号的产生
%y=pulstran(t,d,'func',p1,p2,...)
%其中的p1,p2,...为需要转送给func函数的额外输入参数值（除了变量t之外）

clear all
T=0:1/1E3:1;
D=0:1/4:1;
Y=pulstran(T,D,'rectpuls',0.1);%生成周期性的矩形脉冲信号
subplot(121)
plot(T,Y)
xlabel('t')
ylabel('h(t)')
grid on
axis([0 1 -0.1 1.1])
T=0:1/1E3:1;
D=0:1/3:1;
Y=pulstran(T,D,'tripuls',0.2,1);%生成周期性的三角波信号
subplot(122)
plot(T,Y)
xlabel('t')
ylabel('h(t)')
grid on
axis([0 1 -0.1 1.1])
