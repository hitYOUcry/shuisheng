%AWGN phase unkown
clear
clc

A=1;        %信号幅度
fc=240;     %信号频率
fai=pi/4;   %初始相位
fs=4000;    %采样频率
T=0.01;     %观测时间
N=T*fs;     %数据点数

n=1:N;
S=A*sin(2*pi*fc*n/fs+fai);  %信号
Ps=0.5*A^2;                 %信号功率
SNR=-5;          %信噪比
Pn=Ps/10^(SNR/10);  %噪声功率
N0=2*Pn;

%计算门限
Pf=[0.005:0.005:0.045 0.05:0.05:1];
% Pf=0.8;
Z0=sqrt(-log(Pf)*N0*Ps*N);

Pd_all=[];
for m=1:length(Pf)
    exisit=0;       
    number=5000;    %实验次数
    for i=1:number
        noise=sqrt(Pn)*randn(1,N);
        x=1.*S+noise;        
        
        Gs=sum(x.*A.*sin(2*pi*fc*n/fs));
        Gc=sum(x.*A.*cos(2*pi*fc*n/fs));
        Z=sqrt(Gs^2+Gc^2);
        if Z>Z0(m)
            exisit=exisit+1;
        end    
    end
    Pd=exisit/number;
    Pd_all=[Pd_all Pd];
end
hold on;
plot([0 Pf],[0 Pd_all]);
axis([0 1 0 1.1]);
grid on;

