%AWGN phase unkown
clear
clc

A=1;        %�źŷ���
fc=240;     %�ź�Ƶ��
fai=pi/4;   %��ʼ��λ
fs=4000;    %����Ƶ��
T=0.01;     %�۲�ʱ��
N=T*fs;     %���ݵ���

n=1:N;
S=A*sin(2*pi*fc*n/fs+fai);  %�ź�
Ps=0.5*A^2;                 %�źŹ���
SNR=-5;          %�����
Pn=Ps/10^(SNR/10);  %��������
N0=2*Pn;

%��������
Pf=[0.005:0.005:0.045 0.05:0.05:1];
% Pf=0.8;
Z0=sqrt(-log(Pf)*N0*Ps*N);

Pd_all=[];
for m=1:length(Pf)
    exisit=0;       
    number=5000;    %ʵ�����
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

