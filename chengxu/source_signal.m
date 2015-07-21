%% 信号源仿真 %%
clc;
close all;


L = 50; % 快拍数
M = 8; % 阵元数
k = 2; % 信源数
SNR = 5; % 信噪比
d = 0.5; %
t = 0:L-1;
K = sqrt(2*10^(SNR/10)); % 信号幅度
the0 = [40 45];
the1 = the0(1) * pi /180;
the2 = the0(2) * pi /180;

%% 产生源信0号
w = [pi/8,pi/6]; % frequency
fai1 = unidrnd(360,1,L)*pi/180;
fai2 = unidrnd(360,1,L)*pi/180;

s1 = K*exp(1i * (w(1) * t + fai1));
s2 = K*exp(1i * (w(2) * t + fai2));
figure;
subplot(1,2,1);
plot(real(K*exp(1i * (w(1) * t))));
title('w = pi/8 的正弦信号');
xlabel('采样点');
ylabel('幅值');

subplot(1,2,2);
scatter(t,fai1,'+');
title('随机相位');
xlabel('采样点');
ylabel('相位（pi)')

figure;
plot(real(s1));
title('信源一的波形图');
xlabel('采样点');
ylabel('幅值');

%% 高斯白噪声
Nn = wgn(M,L,0,'complex');
figure;
plot(t,real(Nn(1,:)));
title('高斯白噪声');
xlabel('采样点');
ylabel('幅值');

%% 阵列接收信号
m = 0:M-1;
Aerfa1 = exp(-1i * 2 * pi * (m)' * d * sin(the1));
Aerfa2 = exp(-1i * 2 * pi * (m)' * d * sin(the2));
A = [Aerfa1 Aerfa2];
ss = [s1;s2];
x = A * ss + Nn;
figure;
plot(real(x(1,:)));
title('阵元一');
xlabel('采样点');
ylabel('幅值');
