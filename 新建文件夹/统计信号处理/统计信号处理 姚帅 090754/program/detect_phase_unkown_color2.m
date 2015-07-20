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
T=0.03;     %�۲�ʱ��
fs=4000;
N=T*fs;     %���ݵ���
n=1:N;
S=A*sin(2*pi*fc*n/fs+fai);  %�ź�


nr=0:4;
Rr = 2*exp(-wm*nr/fs).*cos(w0*nr/fs)+2*K*sin(w0*nr/fs).*exp(-wm*nr/fs);

SNR=-5;          %�����
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

Ss=A*sin(2*pi*fc*n/fs);
Sc=A*cos(2*pi*fc*n/fs);

%û���ź�ʱGs��Gc�ķ��������ȣ�ȡƽ��
D_Gs=Ss*C_inv*Ss';
D_Gc=Sc*C_inv*Sc';
D=(D_Gs+D_Gc)/2   %�����ֲ��ķ���

%ȷ������
% Pf=0.1;        %�龯����
Pf=[0.005:0.005:0.045 0.05:0.05:1];
% Pf=[0.005:0.005:0.3];
% Pf=[0.05:0.05:1];  %�龯����
Z0=icdf('rayl',1-Pf,sqrt(D));%�ñ�׼��
% Z0=icdf('rayl',Pf,D);%�÷���

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

        %����ͳ����
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
title('ROC����(��ɫ��˹���������������λ�źŵļ��)');
xlabel('�龯����');
ylabel('������');
axis([0 1 0 1.1]);
grid on;
t1=toc

