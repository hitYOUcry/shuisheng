%��ɫ�����źŵļ�⣬��λƵ��δ֪  ����ķ���
clear
clc

K=0.1;
fm=37400;
wm=2*pi*fm;
f0=37600;
w0=2*pi*f0;

A=1;        %�źŷ���
T=0.0001;     %�۲�ʱ��
fs=1000000;
N=T*fs;     %���ݵ���
n=1:N;

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

%ȷ������
% Pf=0.1;        %�龯����
Pf=[0.005:0.005:0.045 0.05:0.05:1];


Pd_all=[];
tic;
for m=1:length(Pf)
    exisit=0;
    Z=[];
    number=1000;
    fai=2*pi*rand(1,number);
    f=40*rand(1,number)+220;
    f_temp=220:260;
    fL=length(f_temp);    
    
    for i=1:number        
        S=A*sin(2*pi*f(i)*n/fs+fai(i));  %�ź�
        g=sigma*randn(1,N);
        noise=filter(b0,a,g);
        x=1.*S+noise;

        for j=1:fL
            Ss(j,:)=A*sin(2*pi*f_temp(j)*n/fs);
            Sc(j,:)=A*cos(2*pi*f_temp(j)*n/fs);

            %û���ź�ʱGs��Gc�ķ��������ȣ�ȡƽ��
            D_Gs(j)=Ss(j,:)*C_inv*Ss(j,:)';
            D_Gc(j)=Sc(j,:)*C_inv*Sc(j,:)';
            D(j)=(D_Gs(j)+D_Gc(j))/2;   %�����ֲ��ķ���
            Z0(j)=icdf('rayl',1-Pf(m),sqrt(D(j)));%�ñ�׼��
            %����ͳ����
            Gs(j)=x*C_inv*Ss(j,:)';
            Gc(j)=x*C_inv*Sc(j,:)';
            Z(i,j)=sqrt(Gs(j)^2+Gc(j)^2);
            Z_guiyi(i,j)=Z(i,j)/Z0(j);

        end
        
%         [ZZ(i) index(i)]=max(Z_guiyi(i,:));
%         f_match=min(f_temp)+index-1;
        if max(Z_guiyi(i,:))>1
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

