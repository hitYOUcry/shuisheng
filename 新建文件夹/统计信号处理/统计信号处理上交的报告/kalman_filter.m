function []=kalman_filter()
%%%kalman_filter: 实现卡尔曼滤波得目标的运动轨迹
%%%
%%%
%%%by tanjunhong,on 7.6,2011
%%%modified on


clc
clear all
% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% set parameter %%%%%%%%%%
distance0=2e4;          %目标初始距离
doa0=45;                %目标初始方位
v0=20;                  %目标径向速度，大约40节/时，1节=1海里=1.852公里
w0=1/60;                %目标角速度，度/秒
T=2;                    %观测时间间隔
time=2*2600;            %观测时间
number=time/T;

w1=-0.8+(0.8+0.8)*rand(1,number);
w2=-0.001+(0.001+0.001)*rand(1,number);
W=[w1;w2];               %状态噪声矢量矩阵2*number
n1=500*randn(1,number);
n2=2*randn(1,number);
N=[n1;n2];               %观测噪声矢量矩阵

% var_w1=w1*w1'/number
% var_n1=n1*n1'/number
% var_w2=w2*w2'/number
% var_n2=n2*n2'/number
% Q=[var_w1,0;0,var_w2];    
% R=[var_n1,0;0,var_n2];
Q=[0.8^2/3,0;0,0.001^2/3]; %状态噪声(均匀分布)协方差矩阵       ？？？？方差与T有关？？？？
R=[500^2,0;0,2^2];         %观测噪声协方差矩阵

Fai=[1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];
G=[T^2/2 0;T 0;0 T^2/2;0 T];
C=[1 0 0 0;0 0 1 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 无状态噪声和观测噪声x_1=s_1   %%
s_1(:,1)=[distance0;v0;doa0;w0];    %信号模型初始化【s1,s2,s3,s4】s1---k时刻目标径向距离；s2---k时刻目标径向速度
                                    %                            s3---k时刻目标方位；    s4---k时刻目标角速度
for k=2:number
    s_1(:,k)=Fai*s_1(:,k-1);
end
figure
subplot(2,2,1)
polar(s_1(3,:)*pi/180,s_1(1,:),'b');
title('无状态噪声、无观测噪声');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 有状态噪声x_2=s_2， 有观测噪声 s_3 %%%
s_2(:,1)=[distance0;v0;doa0;w0];
for k=2:number
    s_2(:,k)=Fai*s_2(:,k-1)+G*W(:,k);    %信号模型（含状态噪声）
end
x=C*s_2+N;                          %含状态噪声和观测噪声的观测模型【x1,x2】x1--k时刻目标距离测量值，x2--k时刻方位测量值
subplot(2,2,2);
polar(s_2(3,:)*pi/180,s_2(1,:),'b');
title('有状态噪声、无观测噪声'); %有状态噪声、有观测噪声时目标的实际轨迹（排除观测噪声的干扰）
subplot(2,2,3);
polar(x(2,:)*pi/180,x(1,:),'b');
title('有状态噪声、有观测噪声');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%卡尔曼滤波，估计距离和方位 %%%
s_est(:,1)=[distance0;v0;doa0;w0];
P(:,:,1)=s_est(:,1)*s_est(:,1)'/4;
for k=2:number
    
    s_temp=Fai*s_est(:,k-1);              %预测方程
    x_temp=x(:,k)-C*s_temp;               %新息方程
    P_temp=Fai*P(:,:,k-1)*Fai'+G*Q*G';    %协方差方程
    B_temp=P_temp*C'*(R+C*P_temp*C')^(-1); %增益矩阵方程
    P(:,:,k)=(eye(4)-B_temp*C)*P_temp;
    s_est(:,k)=s_temp+B_temp*x_temp;      %滤波方程
    
end

subplot(2,2,4)
polar(s_est(3,:)*pi/180,s_est(1,:),'b');
title('估计轨迹');


figure
subplot(2,2,1)
plot(s_1(3,:),s_1(1,:)/1000,'b');
xlabel('目标方位 \theta (^o)');ylabel('目标径向距离 （km）')
title('无状态噪声、无观测噪声');
subplot(2,2,2);
plot(s_2(3,:),s_2(1,:)/1000,'b');
xlabel('目标方位 \theta (^o)');ylabel('目标径向距离 （km）')
title('有状态噪声、无观测噪声'); %有状态噪声、有观测噪声时目标的实际轨迹（排除观测噪声的干扰）
subplot(2,2,3);
plot(x(2,:),x(1,:)/1000,'b');
xlabel('目标方位 \theta (^o)');ylabel('目标径向距离 （km）')
title('有状态噪声、有观测噪声');
subplot(2,2,4)
plot(s_est(3,:),s_est(1,:)/1000,'b');
xlabel('目标方位 \theta (^o)');ylabel('目标径向距离 （km）')
title('有状态噪声、有观测噪声');


figure
subplot(2,2,1)
plot(s_1(1,:)/1000,s_1(3,:),'b');
ylabel('目标方位 \theta (^o)');xlabel('目标径向距离 （km）')
title('无状态噪声、无观测噪声');
subplot(2,2,2);
plot(s_2(1,:)/1000,s_2(3,:),'b');
ylabel('目标方位 \theta (^o)');xlabel('目标径向距离 （km）')
title('有状态噪声、无观测噪声'); %有状态噪声、有观测噪声时目标的实际轨迹（排除观测噪声的干扰）
subplot(2,2,3);
plot(x(1,:)/1000,x(2,:),'b');
ylabel('目标方位 \theta (^o)');xlabel('目标径向距离 （km）')
title('有状态噪声、有观测噪声');
subplot(2,2,4)
plot(s_est(1,:)/1000,s_est(3,:),'b');
ylabel('目标方位 \theta (^o)');xlabel('目标径向距离 （km）')
title('估计轨迹');


    
    

    


               









