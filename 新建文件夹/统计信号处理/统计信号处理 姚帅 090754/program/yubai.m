%��ɫ�����źŵļ�⣬��λδ֪  ����ķ���
clear
clc

K=0.1;
fm=200;
wm=2*pi*fm;
f0=260;
w0=2*pi*f0;

A=1;        %�źŷ���
fc=240;     %�ź�Ƶ��
fai=pi/4;   %��ʼ��λ
T=0.2;     %�۲�ʱ��
fs=4000;
N=T*fs;     %���ݵ���
n=1:N;
S=A*sin(2*pi*fc*n/fs+fai);  %�ź�


nr=0:4;
Rr = 2*exp(-wm*nr/fs).*cos(w0*nr/fs)+2*K*sin(w0*nr/fs).*exp(-wm*nr/fs);

SNR=0;          %�����
Ps=0.5*A^2;                 %�źŹ���
Pn=Ps/10^(SNR/10);  %��������
sigma2=Pn/Rr(1);
sigma=sqrt(sigma2);

% Rr=Rr/max(Rr);
[a,epsilon]=rtoa(Rr);
b0=sqrt(epsilon);

%���������
r=[];
r(1:5)=Rr;
for k=6:N
    r(k)=-r(k-4:k-1)*fliplr(a(2:5)')';
end
C=toeplitz(r).*sigma2;      %��ɫ������Э������
C_inv=C^(-1);
DD=chol(C_inv);

P_noise=0;
P_bai=0;
times=1000;
for i=1:times
g=sigma*randn(1,N);
noise=filter(b0,a,g);

bai=noise*DD;
P_noise=P_noise+abs(fft(noise)).^2/N;
P_bai=P_bai+abs(fft(bai)).^2/N;
end
f=n*fs/N;
plot(f,P_noise/times);
hold on;
plot(f,P_bai/times,'r');

