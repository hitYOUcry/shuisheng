%AWGN phase frequency unkown
clear
clc

tic;
A=1;        %�źŷ���
% fc=240;     %�ź�Ƶ��
% fai=pi/4;   %��ʼ��λ
fs=4000;    %����Ƶ��
T=0.015;     %�۲�ʱ��
N=T*fs;     %���ݵ���

n=1:N;
% S=A*sin(2*pi*fc*n/fs+fai);  %�ź�
Ps=0.5*A^2;                 %�źŹ���
SNR=-5;          %�����
Pn=Ps/10^(SNR/10);  %��������
N0=2*Pn;

%��������
Pf=[0.005:0.005:0.045 0.05:0.05:1];
% Pf=0.01;
Z0=sqrt(-log(Pf)*N0*Ps*N);

Pd_all=[];
for m=1:29
% for m=1
    exisit=0;       
    number=4000;    %ʵ�����
    fai=2*pi*rand(1,number);
    f=40*rand(1,number)+220;
    
    for i=1:number
%     for i=1
        S=A*sin(2*pi*f(i)*n/fs+fai(i));  %�ź�
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
