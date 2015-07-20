%有色噪声信号的检测，相位未知  矩阵的方法
clear
clc

K=0.1;
fm=200;
wm=2*pi*fm;
f0=260;
w0=2*pi*f0;

A=1;        %信号幅度
fc=240;     %信号频率
fai=pi/4;   %初始相位
T=0.03;     %观测时间
fs=4000;
N=T*fs;     %数据点数
n=1:N;
S=A*sin(2*pi*fc*n/fs+fai);  %信号


nr=0:4;
Rr = 2*exp(-wm*nr/fs).*cos(w0*nr/fs)+2*K*sin(w0*nr/fs).*exp(-wm*nr/fs);

SNR=-5;          %信噪比
Ps=0.5*A^2;                 %信号功率
Pn=Ps/10^(SNR/10);  %噪声功率
sigma2=Pn/Rr(1);
sigma=sqrt(sigma2);

% Rr=Rr/max(Rr);
[a,epsilon]=rtoa(Rr);
b0=sqrt(epsilon);

%自相关外推
r=[];
r(1:5)=Rr;
for k=6:N
    r(k)=-r(k-4:k-1)*fliplr(a(2:5)')';
end
C=toeplitz(r).*sigma2;      %有色噪声的协方差阵
C_inv=C^(-1);
DD=chol(C_inv);

Ss=A*sin(2*pi*fc*n/fs);
Sc=A*cos(2*pi*fc*n/fs);

%没有信号时Gs与Gc的方差，近似相等，取平均
D_Gs=Ss*C_inv*Ss';
D_Gc=Sc*C_inv*Sc';
D=(D_Gs+D_Gc)/2   %瑞利分布的方差

%确定门限
% Pf=0.1;        %虚警概率
Pf=[0.005:0.005:0.045 0.05:0.05:1];
% Pf=[0.005:0.005:0.3];
% Pf=[0.05:0.05:1];  %虚警概率
Z0=icdf('rayl',1-Pf,sqrt(D));%用标准差
% Z0=icdf('rayl',Pf,D);%用方差

Pd_all=[];
tic;
for m=1:length(Pf)
    exisit=0;
    Z=[];

    number=5000;
    for i=1:number
        g=sigma*randn(1,N);
        noise=filter(b0,a,g);
        x=1.*S+noise;

        %计算统计量
        Gs=x*C_inv*Ss';
        Gc=x*C_inv*Sc';
        Z(i)=sqrt(Gs^2+Gc^2);

        if Z(i)>Z0(m)
            exisit=exisit+1;
        end
    end
    Pd=exisit/number;
    Pd_all=[Pd_all Pd];
end
hold on;
plot([0 Pf],[0 Pd_all]);
title('ROC曲线(有色高斯噪声背景下随机相位信号的检测)');
xlabel('虚警概率');
ylabel('检测概率');
axis([0 1 0 1.1]);
grid on;
t1=toc

