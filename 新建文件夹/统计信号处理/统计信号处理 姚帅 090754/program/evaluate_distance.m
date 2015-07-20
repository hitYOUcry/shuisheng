clear
clc

K=0.1;
fm=200;
wm=2*pi*fm;
f0=260;
w0=2*pi*f0;

A=1;        %信号幅度
fc=200;     %信号频率
fai=pi/4;   %初始相位
T_pulse=0.2;     %脉宽
T=5;       %观测时间
fs=4000;
N=T*fs;     %数据点数
N_pulse=T_pulse*fs;         %发出的正弦波的宽度(点数)
T_win=0.02;
N_win=T_win*fs;              %窗口宽度(点数)
n=1:N_pulse;
S=A*sin(2*pi*fc*n/fs+fai);  %脉冲信号

distance=1300; %目标与第一个基元的位置
c=1500;         %声速
tao=2*distance/c;
tao_n=tao*fs;

nr=0:4;
Rr = 2*exp(-wm*nr/fs).*cos(w0*nr/fs)+2*K*sin(w0*nr/fs).*exp(-wm*nr/fs);

SNR=-5;          %信噪比
Ps=0.5*A^2;                 %信号功率
Pn=Ps/10^(SNR/10);  %噪声功率
sigma2=Pn/Rr(1);
sigma=sqrt(sigma2);

[a,epsilon]=rtoa(Rr);
b0=sqrt(epsilon);

%定门限
%自相关外推
r=[];
r(1:5)=Rr;
for k=6:N_win
    r(k)=-r(k-4:k-1)*fliplr(a(2:5)')';
end
C=toeplitz(r).*sigma2;      %噪声的协方差阵
k=1:N_win;
Ss=A*sin(2*pi*fc*k/fs);
Sc=A*cos(2*pi*fc*k/fs);
%没有信号时Gs与Gc的方差，近似相等，取平均
D_Gs=Ss*C^(-1)*Ss';
D_Gc=Sc*C^(-1)*Sc';
D=(D_Gs+D_Gc)/2;    %瑞利分布的方差
%%%%%%    10^(SNR/10)*N_win=Ps*N_win/Pn
Pf=0.01;  %虚警概率
Z0=icdf('rayl',1-Pf,sqrt(D));%用标准差

g=sigma*randn(1,N);
noise=filter(b0,a,g);

x_all=noise;
x_all(round(tao_n)+1:round(tao_n)+N_pulse)=x_all(round(tao_n)+1:round(tao_n)+N_pulse)+S;
plot(x_all)
Number_win=N/N_win;
exisit=zeros(1,Number_win);
for i=1:Number_win
    x=x_all(1+(i-1)*N_win:i*N_win);
    %计算统计量
    Gs=x*C^(-1)*Ss';
    Gc=x*C^(-1)*Sc';
    Z=sqrt(Gs^2+Gc^2);
    if Z>Z0
        exisit(i)=exisit(i)+1;
    end
end
figure;
stem(exisit);
for i=1:Number_win-N_win
    temp=sum(exisit(i:i+N_pulse/N_win-1));
    if (exisit(i)==1)&&(temp>5)
        index=i
        break
    end
end

% exisit2=zeros(1,N_win-1);
% for j=1:N_win-1
%     x=x_all(j+1+(index-2)*N_win:(index-1)*N_win+j);
%     %计算统计量
%     Gs=x*C^(-1)*Ss';
%     Gc=x*C^(-1)*Sc';
%     Z=sqrt(Gs^2+Gc^2);
%     if Z>Z0
%         exisit2(j)=exisit2(j)+1;
%     end    
% end
% for j=1:N_win-10
%     temp=sum(exisit2(j:j+8));
%     if (exisit2(j)==1)&&(temp>5)
%         index2=j
%         break
%     end
% end
% figure;
% stem(exisit2);

distance
d_eva=((index-1)*T_win)/2*c
% d_eva_ok=((index-1)*T_win-(length(exisit2)-index2)/fs)/2*c




