%求AR模型的系数,产生混有有色噪声的信号
clear
clc

K=0.1;
fm=200;
wm=2*pi*fm;
f0=260;
w0=2*pi*f0;

%自己算的
fs=4000;
nr=0:4;
Rr = 2*exp(-wm*nr/fs).*cos(w0*nr/fs)+2*K*sin(w0*nr/fs).*exp(-wm*nr/fs);
Rr=Rr/max(Rr);
[a,epsilon]=rtoa(Rr);
[H,TT]=impz(1,a);


A=1;        %信号幅度
fc=240;     %信号频率
fai=pi/4;   %初始相位
T=1;     %观测时间
N=T*fs;     %数据点数
n=1:N;
S=A*sin(2*pi*fc*n/fs+fai);  %信号
Ps=0.5*A^2;                 %信号功率
SNR=-15;          %信噪比
Pn=Ps/10^(SNR/10);  %噪声功率
sigma=Pn/sum(H.^2);

g=sqrt(sigma)*randn(1,N);
noise=filter(1,a,g);

x=1.*S+noise;
% plot(abs(fft(x)).^2/N)

plot([0:3999],abs(fft(x)))
axis([0 2000 0 2000]);

