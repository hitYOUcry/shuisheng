%% 卡尔曼滤波
%--卡尔曼滤波程序
close all;
clc;
x10=50;%目标初始距离20km
x20=0.02;%目标径向速度 20m/s
x30=120;%初始角度
x40=6/60;%角速度  一分钟6°，一小时360°。
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

subplot(2,2,1);
plot(u1);
title('u1: 随机径向加速度');
subplot(2,2,2);
plot(u2);
title('u2: 随机角加速度');
subplot(2,2,3);
plot(v1);
title('v1: 径向测量误差');
subplot(2,2,4);
plot(v2);
title('v2: 方位测量误差');

%产生无状态噪声和观测噪声时由距离和方位表示的目标运动轨迹序列；
x1=zeros(1,N);x2=zeros(1,N);x3=zeros(1,N);x4=zeros(1,N);
x1(1)=x10;
x2(1)=x20;
x3(1)=x30;
x4(1)=x40;
A=[ 1,T,0,0;
    0,1,0,0;
    0,0,1,T;
    0,0,0,1];
C=[ 1,0,0,0;
    0,0,1,0];
for i=1:(N-1),
    x_pre=[x1(i);x2(i);x3(i);x4(i)];
    x=A*x_pre;
    x1(i+1)=x(1);
    x2(i+1)=x(2);
    x3(i+1)=x(3);
    x4(i+1)=x(4);
end
z1=x1;
z2=x3;

figure;
subplot(2,2,1);
polar((z2/180)*pi,z1);
title('无观测和状态噪声时运动轨迹');
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
subplot(2,2,2);
polar((z2/180)*pi,z1);
title('有状态噪声时运动轨迹');
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
subplot(2,2,3);
polar((z2/180)*pi,z1);
title('有观测和状态噪声时运动轨迹');

%卡尔曼滤波
R=[p_v1^2,0;0,p_v2^2];
Q=[(T^4/4)*p_u1_2,  (T^3/2)*p_u1_2,             0,              0;
   (T^3/2)*p_u1_2,    (T^2)*p_u1_2,             0,              0;
    0           ,        0,       (T^4/4)*p_u2_2,   (T^3/2)*p_u2_2;
    0           ,        0,       (T^3/2)*p_u2_2,     (T^2)*p_u2_2];
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
subplot(2,2,4);
polar((x(3,:)/180)*pi,x(1,:));    
title('卡尔曼估计');
angle2=x(3,:);
distance2=x(1,:);
figure;
polar((angle1/180)*pi,distance1,'b');
hold on;
polar((angle2/180)*pi,distance2,'r');
title('原轨迹和卡尔曼滤波轨迹比较');
legend('原轨迹','卡尔曼估计',4);
hold off;