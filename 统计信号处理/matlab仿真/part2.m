%% 卡尔曼滤波 目标轨迹跟踪

close all;
clear;
clc;

%% 值初始化
target_dis = 50;% km
radial_velocity = 0.02;% m/s
initAngle = 120;
angular_velocity = 6/60;
T = 2;
N = 1800;

p_u1 = 0.0002;%2m/s2
p_u2 = 0.01/60;%一分钟0.1°的加速度
%观测噪声方差
p_v1 = 5;
p_v2 = 4;

%% 根据噪声模型产生状态噪声和观测噪声序列
%径向距离加速度
u1 = p_u1-2*p_u1*rand(1,N);
p_u1_2 = p_u1^2/3;
u2 = p_u2-2*p_u2*rand(1,N);%角加速度
p_u2_2 = p_u2^2/3;
v1 = p_v1*randn(1,N);%距离测量误差
v2 = p_v2*randn(1,N);%方位测量误差


%% 运动目标轨迹；

A = [ 1,T,0,0;
    0,1,0,0;
    0,0,1,T;
    0,0,0,1];
C = [ 1,0,0,0;
    0,0,1,0];
D = [(T^2)/2,0;
    T,0;
    0,(T^2)/2;
    0,T];
x1 = zeros(1,N);x2 = zeros(1,N);x3 = zeros(1,N);x4 = zeros(1,N);
x1(1) = target_dis;
x2(1) = radial_velocity;
x3(1) = initAngle;
x4(1) = angular_velocity;
for i = 1:(N-1),
    x_pre2 = [x1(i);x2(i);x3(i);x4(i)];
    U = [u1(i);u2(i)];
    xx2 = A*x_pre2+D*U;
    x1(i+1) = xx2(1);
    x2(i+1) = xx2(2);
    x3(i+1) = xx2(3);
    x4(i+1) = xx2(4);
end
z1 = x1;
z2 = x3;
figure;
polar((z2/180)*pi,z1);
title('运动目标轨迹');

%% 运动目标轨迹(含噪)

D = [(T^2)/2,0;
    T,0;
    0,(T^2)/2;
    0,T];
x1 = zeros(1,N);x2 = zeros(1,N);x3 = zeros(1,N);x4 = zeros(1,N);
x1(1) = target_dis;
x2(1) = radial_velocity;
x3(1) = initAngle;
x4(1) = angular_velocity;
for i = 1:(N-1),
    x_pre2 = [x1(i);x2(i);x3(i);x4(i)];
    U = [u1(i);u2(i)];
    xx2 = A*x_pre2+D*U;
    x1(i+1) = xx2(1);
    x2(i+1) = xx2(2);
    x3(i+1) = xx2(3);
    x4(i+1) = xx2(4);
end
z1 = x1;
z2 = x3;
angle1 = x3;
distance1 = x1;
%根据观测噪声和观测方程产生观测序列，观测为距离和方位；
for i = 1:(N-1),
    x_pre = [x1(i);x2(i);x3(i);x4(i)];
    U = [u1(i);u2(i)];
    x = A*x_pre+D*U;
    x1(i+1) = x(1);
    x2(i+1) = x(2);
    x3(i+1) = x(3);
    x4(i+1) = x(4);
end
z1 = x1+v1;
z2 = x3+v2;
figure;
polar((z2/180)*pi,z1);
title('含噪的运动目标轨迹');

%% 卡尔曼滤波
R = [p_v1^2,0;0,p_v2^2];
Q = [(T^4/4)*p_u1_2,(T^3/2)*p_u1_2,0,0;
   (T^3/2)*p_u1_2,(T^2)*p_u1_2,0,0;
    0,0,(T^4/4)*p_u2_2,(T^3/2)*p_u2_2;
    0,0,(T^3/2)*p_u2_2,(T^2)*p_u2_2];
x(:,1) = [x1(1);x2(1);x3(1);x4(1)];
P = x(1)*(x(1))';
for k = 1:N,
  y(:,k) = [z1(k);z2(k)];
end;
I = [1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,1];
for k = 1:(N-1),
  xk_1 = A*x(:,k);   
  p_tran = A*P*A'+Q;
  K1(:,:,k) = p_tran*C'*inv(C*p_tran*C'+R);%卡尔曼滤波系数K
  x(:,k+1) = xk_1+K1(:,:,k)*(y(:,k+1)-C*xk_1);
  P = (I-K1(:,:,k)*C)*p_tran;
end
figure;
polar((x(3,:)/180)*pi,x(1,:));    
title('卡尔曼估计');
angle2 = x(3,:);
distance2 = x(1,:);
figure;
polar((angle1/180)*pi,distance1,'b');
hold on;
polar((angle2/180)*pi,distance2,'r');
title('原轨迹和卡尔曼滤波跟踪轨迹比较');
legend('原轨迹','卡尔曼估计',4);
hold off;