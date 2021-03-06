%高精度延时滤波器
clear
clc

A=1;
f=200;
fs=4000;
fai=pi/4;
T=0.5;
N=fs*T;
n=1:N;
s=A*sin(2*pi*f*n/fs);

M=12;
k=0:M;
taof=0.5;       %taof的取值范围是-0.5~0.5
hd=sin(k*pi-taof*pi)./(k*pi-taof*pi);
% win=ones(1,M+1);
win=hanning(M+1)';  %矩形窗效果不好，用hanning窗

g=hd.*win;
h=g/sum(g);

%用正弦信号做高精度延时
s_after=conv(s,h);

figure
plot(s,'r.');    %s(t)
hold on;
plot(s_after,'b.'); %s(t-taof)
s_true=A*sin(2*pi*f*(n-taof)/fs);
hold on;
plot(s_true,'g*');%延时后的真值


% %用白噪声做高精度延时
% s1=A*sin(2*pi*f*n/fs+fai);
% % s1=randn(1,N);
% for i=1:N/2    
%     s2(i)=s1(i*2);        %偶数点
%     s3(i)=s1(i*2-1);      %奇数点
% end
% safter=conv(s3,h);        %延迟0.5个采样周期
% figure;
% plot(s2,'r');
% hold on;
% plot(safter);
% 



% %有色噪声
% K=0.1;
% fm=200;
% wm=2*pi*fm;
% f0=260;
% w0=2*pi*f0;
% nr=0:4;
% Rr = 2*exp(-wm*nr/fs).*cos(w0*nr/fs)+2*K*sin(w0*nr/fs).*exp(-wm*nr/fs);
% [a,epsilon]=rtoa(Rr);
% b0=sqrt(epsilon);
% 
% SNR=-5;          %信噪比
% Ps=0.5*A^2;                 %信号功率
% Pn=Ps/10^(SNR/10);  %噪声功率
% sigma2=Pn/Rr(1);
% sigma=sqrt(sigma2);
% 
% g=sigma*randn(1,N);
% noise=filter(b0,a,g);
% s1=1.*s+noise;
% 
% for i=1:N/2    
%     s2(i)=s1(i*2);        %偶数点
%     s3(i)=s1(i*2-1);      %奇数点
% end
% safter=conv(s3,h);        %延迟0.5个采样周期
% figure;
% plot(s3,'r');
% hold on;
% plot(s2,'g');
% hold on;
% plot(safter);

