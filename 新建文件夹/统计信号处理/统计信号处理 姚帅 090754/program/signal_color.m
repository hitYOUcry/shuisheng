%��ARģ�͵�ϵ��,����������ɫ�������ź�
clear
clc

K=0.1;
fm=200;
wm=2*pi*fm;
f0=260;
w0=2*pi*f0;

%�Լ����
fs=4000;
nr=0:4;
Rr = 2*exp(-wm*nr/fs).*cos(w0*nr/fs)+2*K*sin(w0*nr/fs).*exp(-wm*nr/fs);
Rr=Rr/max(Rr);
[a,epsilon]=rtoa(Rr);
[H,TT]=impz(1,a);


A=1;        %�źŷ���
fc=240;     %�ź�Ƶ��
fai=pi/4;   %��ʼ��λ
T=1;     %�۲�ʱ��
N=T*fs;     %���ݵ���
n=1:N;
S=A*sin(2*pi*fc*n/fs+fai);  %�ź�
Ps=0.5*A^2;                 %�źŹ���
SNR=-15;          %�����
Pn=Ps/10^(SNR/10);  %��������
sigma=Pn/sum(H.^2);

g=sqrt(sigma)*randn(1,N);
noise=filter(1,a,g);

x=1.*S+noise;
% plot(abs(fft(x)).^2/N)

plot([0:3999],abs(fft(x)))
axis([0 2000 0 2000]);

