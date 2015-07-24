%% ͳ���źŴ��� ���������ʵ�����


close all;
clear;
clc;

%% ��������
N = 50; 
% ��Ԫ��Ŀ
M = 12; 
%�ź�Դ��Ŀ
k = 2;
% �����
SNR = 5; 
% �źŷ���
Amp = sqrt(2*10^(SNR/10));
%�����
IncAngle = [25 55];
IncAngle1 = IncAngle(1) * pi /180;
IncAngle2 = IncAngle(2) * pi /180;
%�ź�ԴƵ��
w = [pi/6,pi/3];

%% ����Դ�ź�
phase1 = unidrnd(360,1,N)*pi/180;
phase2 = unidrnd(360,1,N)*pi/180;
t = 0:N-1;
s1 = Amp*exp(1i * (w(1) * t + phase1));
s2 = Amp*exp(1i * (w(2) * t + phase2));
figure;plot(real(Amp*exp(1i * (w(1) * t))));
title('�����ź�(w=pi/6,50������)');
xlabel('����');
ylabel('����');
figure;scatter(t,phase1,'*');
title('�����λ');
xlabel('����');
ylabel('��λֵ/pi')
figure;plot(real(s1));
title('�����λ�����źţ���Դһ��');
xlabel('����');
ylabel('����');

%% ����������
noiseSignal = wgn(M,N,0,'complex');
figure;plot(t,real(noiseSignal(1,:)));
title('��˹�������ź�');
xlabel('����');
ylabel('����');

%% ��Ԫ�����źŷ���
d = 0.5; 
m = 0:M-1;
Aerfa1 = exp(-1i * 2 * pi * (m)' * d * sin(IncAngle1));
Aerfa2 = exp(-1i * 2 * pi * (m)' * d * sin(IncAngle2));
A = [Aerfa1 Aerfa2];
ss = [s1;s2];
x = A * ss + noiseSignal;
figure;plot(real(x(1,:)));
title('1����Ԫ�����ź�');
xlabel('����');
ylabel('����');

%% ��λ�ǹ���
% CBF����
R3 = x * x' / N;
theta1 = 0:0.1:90;
N = length(theta1);
for i1 = 1:N
    a_theta1 = exp(-1i * 2 * pi * d * m' * sin(pi * theta1(i1)/180));
    P_cbf3(i1) = a_theta1' * R3 * a_theta1 / (a_theta1' * a_theta1);
end
figure;
plot(theta1,10 * log10(abs(P_cbf3)));
title(strcat('���沨���γ�','(��Ԫ��Ŀ ��',num2str(M),'��)'));
xlabel('�����/deg');
ylabel('�ռ䷽λ��/dB');
grid on;
% MUSIC 
[p3, r3] = eig(R3);
mr3 = fliplr(r3);
mr3 = fliplr(mr3);
mp3 = fliplr(p3);
Us3(:,1:2) = mp3(:,1:2);
Un3(:,:) = mp3(:,3:M);
theta3 = 0:0.1:90;
N = length(theta3);
for i3 = 1:N
    a_theta3 = exp(-1i * 2 * pi * d * m' * sin(theta3(i3) / 180 * pi));
    P_music3(i3) = abs((a_theta3' * a_theta3)/(a_theta3' * Un3 * Un3' * a_theta3));
end
figure;
plot(theta3,10 * log10(P_music3));
title(strcat('MUSIC����','(��Ԫ��Ŀ ��',num2str(M),'��)'));
xlabel('�����/deg');
ylabel('�ռ䷽λ��/dB');

%% �˶�Ŀ����� %%


%% ��������
x10=35;%Ŀ���ʼ����50km
x20=0.02;%Ŀ�꾶���ٶ� 20m/s
x30=45;
x40=6/60;%���ٶ�
T=2;
N=1800;
%״̬������Χ
p_u1=0.0002;%2m/s2
p_u2=0.01/60;%һ����0.1��ļ��ٶ�
%�۲���������
p_v1=5;
p_v2=4;

%��������ģ�Ͳ���״̬�����͹۲���������
%���������ٶ�
u1=p_u1-2*p_u1*rand(1,N);
p_u1_2=p_u1^2/3;
%�Ǽ��ٶ�
u2=p_u2-2*p_u2*rand(1,N);
p_u2_2=p_u2^2/3;
%����������
v1=p_v1*randn(1,N);
%��λ�������
v2=p_v2*randn(1,N);

%% �˶�Ŀ��켣ģ��
A=[ 1,T,0,0;
    0,1,0,0;
    0,0,1,T;
    0,0,0,1];
C=[ 1,0,0,0;
    0,0,1,0];
D=[(T^2)/2,0;
    T,0;
    0,(T^2)/2;
    0,T];
x1=zeros(1,N);x2=zeros(1,N);x3=zeros(1,N);x4=zeros(1,N);
x1(1)=x10;
x2(1)=x20;
x3(1)=x30;
x4(1)=x40;
for i=1:(N-1),
    x_pre2=[x1(i);x2(i);x3(i);x4(i)];
    U=[u1(i);u2(i)];
    xx2=A*x_pre2+D*U;
    x1(i+1)=xx2(1);
    x2(i+1)=xx2(2);
    x3(i+1)=xx2(3);
    x4(i+1)=xx2(4);
end
z1=x1;
z2=x3;
figure;
polar((z2/180)*pi,z1);
title('�˶�Ŀ��켣');

%% ģ�⺬����˶�Ŀ��켣
%��״̬����ʱ��Ŀ���˶��켣���У�
D=[(T^2)/2,0;
    T,0;
    0,(T^2)/2;
    0,T];
x1=zeros(1,N);x2=zeros(1,N);x3=zeros(1,N);x4=zeros(1,N);
x1(1)=x10;
x2(1)=x20;
x3(1)=x30;
x4(1)=x40;
for i=1:(N-1),
    x_pre2=[x1(i);x2(i);x3(i);x4(i)];
    U=[u1(i);u2(i)];
    xx2=A*x_pre2+D*U;
    x1(i+1)=xx2(1);
    x2(i+1)=xx2(2);
    x3(i+1)=xx2(3);
    x4(i+1)=xx2(4);
end
z1=x1;
z2=x3;
angle1=x3;
distance1=x1;
%���ݹ۲������͹۲ⷽ�̲����۲����У��۲�Ϊ����ͷ�λ��
for i=1:(N-1),
    x_pre=[x1(i);x2(i);x3(i);x4(i)];
    U=[u1(i);u2(i)];
    x=A*x_pre+D*U;
    x1(i+1)=x(1);
    x2(i+1)=x(2);
    x3(i+1)=x(3);
    x4(i+1)=x(4);
end
z1=x1+v1;
z2=x3+v2;
figure;
polar((z2/180)*pi,z1);
title('�������ŵ��˶�Ŀ��켣');

%% �������˲�����
R=[p_v1^2,0;0,p_v2^2];
Q=[(T^4/4)*p_u1_2,(T^3/2)*p_u1_2,0,0;
   (T^3/2)*p_u1_2,(T^2)*p_u1_2,0,0;
    0,0,(T^4/4)*p_u2_2,(T^3/2)*p_u2_2;
    0,0,(T^3/2)*p_u2_2,(T^2)*p_u2_2];
x(:,1)=[x1(1);x2(1);x3(1);x4(1)];
P=x(1)*(x(1))';
for k=1:N,
  y(:,k)=[z1(k);z2(k)];
end;
I=[1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,1];
for k=1:(N-1),
  xk_1=A*x(:,k);   
  p_tran=A*P*A'+Q;
  K1(:,:,k)=p_tran*C'*inv(C*p_tran*C'+R);%�������˲�ϵ��K
  x(:,k+1)=xk_1+K1(:,:,k)*(y(:,k+1)-C*xk_1);
  P=(I-K1(:,:,k)*C)*p_tran;
end
figure;
polar((x(3,:)/180)*pi,x(1,:));    
title('�������˲����ٽ��');
angle2=x(3,:);
distance2=x(1,:);
figure;
polar((angle1/180)*pi,distance1,'b');
hold on;
polar((angle2/180)*pi,distance2,'r');
title('��ʵ�켣�뿨�����˲����ٹ켣�Ա�');
legend('��ʵ�켣','�������˲�����',3);
hold off;