clear;
clc;
%目标初始值
s0=20000;       %目标距离(m)
v0=20;          %目标径向速度(m/s)
angle0=30;      %目标初始角度(度)
w0=1/60;        %目标角速度(1度/s)

T=2;            %两次测量的时间间隔
number=1800;    %测量次数

Fai=[1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];
G=[T^2/2 0;T 0;0 T^2/2;0 T];
C=[1 0 0 0;0 0 1 0];

% u1=(2*rand(1,number)-1)*0.8;
% u2=(2*rand(1,number)-1)*0.001;
%离散的u1与u2是服从正态分布
sigma1=sqrt(1.6^2/12);      %均匀分布的标准差
sigma2=sqrt(0.002^2/12);
sigmau1=sigma1/sqrt(T);
sigmau2=sigma2/sqrt(T);
%状态噪声
u1=sigmau1*randn(1,number);
u2=sigmau2*randn(1,number);

%观测噪声
v1=500*randn(1,number);
v2=2*randn(1,number);
u=[u1;u2];
v=[v1;v2];

%不含噪声的状态方程与观测方程
x(:,1)=Fai*[s0 v0 angle0 w0]';
z(:,1)=C*x(:,1);
for k=2:number
    x(:,k)=Fai*x(:,k-1);
    z(:,k)=C*x(:,k);
end

polar(z(2,:)*pi/180,z(1,:),'k');    %观测运动轨迹
figure;
%含噪声的状态方程与观测方程
for k=2:number
    x(:,k)=Fai*x(:,k-1)+G*u(:,k-1);
    z(:,k)=C*x(:,k)+v(:,k);
end

polar(x(3,:)*pi/180,x(1,:),'g');   %真实运动轨迹
figure;
polar(z(2,:)*pi/180,z(1,:),'y');    %观测运动轨迹

%kalma滤波
%噪声的协方差阵为已知条件
Q=[sigmau1^2 0;0 sigmau2^2];    %状态噪声u的协方差阵
R=[500 0;0 2];                  %观测噪声v的协方差阵
%初始化
x0=[z(1,2) 0 z(2,2) 0]';    %可能为z(1,1)与z(2,1)
x_esti(:,1)=x0;
P(:,:,1)=x0*x0';

for k=2:number
    xtemp=Fai*x_esti(:,k-1);            %预测估计方程
    xinxi=z(:,k)-C*xtemp;               %新息
    Ptemp=Fai*P(:,:,k-1)*Fai'+G*Q*G';   
    B=Ptemp*C'*(R+C*Ptemp*C')^(-1);     %增益矩阵方程
    P(:,:,k)=(eye(4)-B*C)*Ptemp;        %协方差方程
    x_esti(:,k)=xtemp+B*xinxi;          %滤波方程
end

figure;
polar(x(3,:)*pi/180,x(1,:),'r-.');      %kalman波形估计结果


