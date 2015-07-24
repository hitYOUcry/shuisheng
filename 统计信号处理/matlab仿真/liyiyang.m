%% 统计信号处理 线列阵仿真实验程序


close all;
clear;
clc;

%% 变量定义
N = 50; 
% 阵元数目
M = 12; 
%信号源数目
k = 2;
% 信噪比
SNR = 5; 
% 信号幅度
Amp = sqrt(2*10^(SNR/10));
%入射角
IncAngle = [25 55];
IncAngle1 = IncAngle(1) * pi /180;
IncAngle2 = IncAngle(2) * pi /180;
%信号源频率
w = [pi/6,pi/3];

%% 产生源信号
phase1 = unidrnd(360,1,N)*pi/180;
phase2 = unidrnd(360,1,N)*pi/180;
t = 0:N-1;
s1 = Amp*exp(1i * (w(1) * t + phase1));
s2 = Amp*exp(1i * (w(2) * t + phase2));
figure;plot(real(Amp*exp(1i * (w(1) * t))));
title('正弦信号(w=pi/6,50个样点)');
xlabel('样点');
ylabel('幅度');
figure;scatter(t,phase1,'*');
title('随机相位');
xlabel('样点');
ylabel('相位值/pi')
figure;plot(real(s1));
title('随机相位正弦信号（信源一）');
xlabel('样点');
ylabel('幅度');

%% 产生白噪声
noiseSignal = wgn(M,N,0,'complex');
figure;plot(t,real(noiseSignal(1,:)));
title('高斯白噪声信号');
xlabel('样点');
ylabel('幅度');

%% 阵元接收信号仿真
d = 0.5; 
m = 0:M-1;
Aerfa1 = exp(-1i * 2 * pi * (m)' * d * sin(IncAngle1));
Aerfa2 = exp(-1i * 2 * pi * (m)' * d * sin(IncAngle2));
A = [Aerfa1 Aerfa2];
ss = [s1;s2];
x = A * ss + noiseSignal;
figure;plot(real(x(1,:)));
title('1号阵元接收信号');
xlabel('样点');
ylabel('幅度');

%% 方位角估计
% CBF方法
R3 = x * x' / N;
theta1 = 0:0.1:90;
N = length(theta1);
for i1 = 1:N
    a_theta1 = exp(-1i * 2 * pi * d * m' * sin(pi * theta1(i1)/180));
    P_cbf3(i1) = a_theta1' * R3 * a_theta1 / (a_theta1' * a_theta1);
end
figure;
plot(theta1,10 * log10(abs(P_cbf3)));
title(strcat('常规波束形成','(阵元数目 ：',num2str(M),'个)'));
xlabel('入射角/deg');
ylabel('空间方位谱/dB');
grid on;
% MUSIC 
[p3, r3] = eig(R3);
mr3 = fliplr(r3);
mr3 = fliplr(mr3);
mp3 = fliplr(p3);
Us3(:,1:2) = mp3(:,1:2);
Un3(:,:) = mp3(:,3:M);
theta3 = 0:0.1:90;
N = length(theta3);
for i3 = 1:N
    a_theta3 = exp(-1i * 2 * pi * d * m' * sin(theta3(i3) / 180 * pi));
    P_music3(i3) = abs((a_theta3' * a_theta3)/(a_theta3' * Un3 * Un3' * a_theta3));
end
figure;
plot(theta3,10 * log10(P_music3));
title(strcat('MUSIC方法','(阵元数目 ：',num2str(M),'个)'));
xlabel('入射角/deg');
ylabel('空间方位谱/dB');

%% 运动目标跟踪 %%


%% 变量定义
x10=35;%目标初始距离50km
x20=0.02;%目标径向速度 20m/s
x30=45;
x40=6/60;%角速度
T=2;
N=1800;
%状态噪声范围
p_u1=0.0002;%2m/s2
p_u2=0.01/60;%一分钟0.1°的加速度
%观测噪声方差
p_v1=5;
p_v2=4;

%根据噪声模型产生状态噪声和观测噪声序列
%径向距离加速度
u1=p_u1-2*p_u1*rand(1,N);
p_u1_2=p_u1^2/3;
%角加速度
u2=p_u2-2*p_u2*rand(1,N);
p_u2_2=p_u2^2/3;
%距离测量误差
v1=p_v1*randn(1,N);
%方位测量误差
v2=p_v2*randn(1,N);

%% 运动目标轨迹模拟
A=[ 1,T,0,0;
    0,1,0,0;
    0,0,1,T;
    0,0,0,1];
C=[ 1,0,0,0;
    0,0,1,0];
D=[(T^2)/2,0;
    T,0;
    0,(T^2)/2;
    0,T];
x1=zeros(1,N);x2=zeros(1,N);x3=zeros(1,N);x4=zeros(1,N);
x1(1)=x10;
x2(1)=x20;
x3(1)=x30;
x4(1)=x40;
for i=1:(N-1),
    x_pre2=[x1(i);x2(i);x3(i);x4(i)];
    U=[u1(i);u2(i)];
    xx2=A*x_pre2+D*U;
    x1(i+1)=xx2(1);
    x2(i+1)=xx2(2);
    x3(i+1)=xx2(3);
    x4(i+1)=xx2(4);
end
z1=x1;
z2=x3;
figure;
polar((z2/180)*pi,z1);
title('运动目标轨迹');

%% 模拟含噪的运动目标轨迹
%有状态噪声时的目标运动轨迹序列；
D=[(T^2)/2,0;
    T,0;
    0,(T^2)/2;
    0,T];
x1=zeros(1,N);x2=zeros(1,N);x3=zeros(1,N);x4=zeros(1,N);
x1(1)=x10;
x2(1)=x20;
x3(1)=x30;
x4(1)=x40;
for i=1:(N-1),
    x_pre2=[x1(i);x2(i);x3(i);x4(i)];
    U=[u1(i);u2(i)];
    xx2=A*x_pre2+D*U;
    x1(i+1)=xx2(1);
    x2(i+1)=xx2(2);
    x3(i+1)=xx2(3);
    x4(i+1)=xx2(4);
end
z1=x1;
z2=x3;
angle1=x3;
distance1=x1;
%根据观测噪声和观测方程产生观测序列，观测为距离和方位；
for i=1:(N-1),
    x_pre=[x1(i);x2(i);x3(i);x4(i)];
    U=[u1(i);u2(i)];
    x=A*x_pre+D*U;
    x1(i+1)=x(1);
    x2(i+1)=x(2);
    x3(i+1)=x(3);
    x4(i+1)=x(4);
end
z1=x1+v1;
z2=x3+v2;
figure;
polar((z2/180)*pi,z1);
title('噪声干扰的运动目标轨迹');

%% 卡尔曼滤波跟踪
R=[p_v1^2,0;0,p_v2^2];
Q=[(T^4/4)*p_u1_2,(T^3/2)*p_u1_2,0,0;
   (T^3/2)*p_u1_2,(T^2)*p_u1_2,0,0;
    0,0,(T^4/4)*p_u2_2,(T^3/2)*p_u2_2;
    0,0,(T^3/2)*p_u2_2,(T^2)*p_u2_2];
x(:,1)=[x1(1);x2(1);x3(1);x4(1)];
P=x(1)*(x(1))';
for k=1:N,
  y(:,k)=[z1(k);z2(k)];
end;
I=[1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,1];
for k=1:(N-1),
  xk_1=A*x(:,k);   
  p_tran=A*P*A'+Q;
  K1(:,:,k)=p_tran*C'*inv(C*p_tran*C'+R);%卡尔曼滤波系数K
  x(:,k+1)=xk_1+K1(:,:,k)*(y(:,k+1)-C*xk_1);
  P=(I-K1(:,:,k)*C)*p_tran;
end
figure;
polar((x(3,:)/180)*pi,x(1,:));    
title('卡尔曼滤波跟踪结果');
angle2=x(3,:);
distance2=x(1,:);
figure;
polar((angle1/180)*pi,distance1,'b');
hold on;
polar((angle2/180)*pi,distance2,'r');
title('真实轨迹与卡尔曼滤波跟踪轨迹对比');
legend('真实轨迹','卡尔曼滤波跟踪',3);
hold off;