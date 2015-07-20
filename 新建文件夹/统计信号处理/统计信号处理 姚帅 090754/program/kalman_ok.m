clear;
clc;
%Ŀ���ʼֵ
s0=20000;       %Ŀ�����(m)
v0=20;          %Ŀ�꾶���ٶ�(m/s)
angle0=30;      %Ŀ���ʼ�Ƕ�(��)
w0=1/60;        %Ŀ����ٶ�(1��/s)

T=2;            %���β�����ʱ����
number=1800;    %��������

Fai=[1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];
G=[T^2/2 0;T 0;0 T^2/2;0 T];
C=[1 0 0 0;0 0 1 0];

% u1=(2*rand(1,number)-1)*0.8;
% u2=(2*rand(1,number)-1)*0.001;
%��ɢ��u1��u2�Ƿ�����̬�ֲ�
sigma1=sqrt(1.6^2/12);      %���ȷֲ��ı�׼��
sigma2=sqrt(0.002^2/12);
sigmau1=sigma1/sqrt(T);
sigmau2=sigma2/sqrt(T);
%״̬����
u1=sigmau1*randn(1,number);
u2=sigmau2*randn(1,number);

%�۲�����
v1=500*randn(1,number);
v2=2*randn(1,number);
u=[u1;u2];
v=[v1;v2];

%����������״̬������۲ⷽ��
x(:,1)=Fai*[s0 v0 angle0 w0]';
z(:,1)=C*x(:,1);
for k=2:number
    x(:,k)=Fai*x(:,k-1);
    z(:,k)=C*x(:,k);
end

polar(z(2,:)*pi/180,z(1,:),'k');    %�۲��˶��켣
figure;
%��������״̬������۲ⷽ��
for k=2:number
    x(:,k)=Fai*x(:,k-1)+G*u(:,k-1);
    z(:,k)=C*x(:,k)+v(:,k);
end

polar(x(3,:)*pi/180,x(1,:),'g');   %��ʵ�˶��켣
figure;
polar(z(2,:)*pi/180,z(1,:),'y');    %�۲��˶��켣

%kalma�˲�
%������Э������Ϊ��֪����
Q=[sigmau1^2 0;0 sigmau2^2];    %״̬����u��Э������
R=[500 0;0 2];                  %�۲�����v��Э������
%��ʼ��
x0=[z(1,2) 0 z(2,2) 0]';    %����Ϊz(1,1)��z(2,1)
x_esti(:,1)=x0;
P(:,:,1)=x0*x0';

for k=2:number
    xtemp=Fai*x_esti(:,k-1);            %Ԥ����Ʒ���
    xinxi=z(:,k)-C*xtemp;               %��Ϣ
    Ptemp=Fai*P(:,:,k-1)*Fai'+G*Q*G';   
    B=Ptemp*C'*(R+C*Ptemp*C')^(-1);     %������󷽳�
    P(:,:,k)=(eye(4)-B*C)*Ptemp;        %Э�����
    x_esti(:,k)=xtemp+B*xinxi;          %�˲�����
end

figure;
polar(x(3,:)*pi/180,x(1,:),'r-.');      %kalman���ι��ƽ��


