%AWGN phase frequency unkown
clear
clc

tic;
A=1;        %信号幅度
% fc=240;     %信号频率
% fai=pi/4;   %初始相位
fs=4000;    %采样频率
T=0.015;     %观测时间
N=T*fs;     %数据点数

n=1:N;
% S=A*sin(2*pi*fc*n/fs+fai);  %信号
Ps=0.5*A^2;                 %信号功率
SNR=-5;          %信噪比
Pn=Ps/10^(SNR/10);  %噪声功率
N0=2*Pn;

%计算门限
Pf=[0.005:0.005:0.045 0.05:0.05:1];
% Pf=0.01;
Z0=sqrt(-log(Pf)*N0*Ps*N);

Pd_all=[];
for m=1:29
% for m=1
    exisit=0;       
    number=4000;    %实验次数
    fai=2*pi*rand(1,number);
    f=40*rand(1,number)+220;
    
    for i=1:number
%     for i=1
        S=A*sin(2*pi*f(i)*n/fs+fai(i));  %信号
        noise=sqrt(Pn)*randn(1,N);
        x=1.*S+noise;
        n=1:N; 
        f_temp=220:260;
        fL=length(f_temp);
        for j=1:fL
            Gs(j)=sum(x.*A.*sin(2*pi*f_temp(j)*n/fs));
            Gc(j)=sum(x.*A.*cos(2*pi*f_temp(j)*n/fs));
        end
        Z=sqrt(Gs.^2+Gc.^2);
        [Zmax index]=max(Z);
        f_match=min(f_temp+index-1);
        if Zmax>Z0(m)
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
t=toc
