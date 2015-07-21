%% �ź�Դ���� %%
clc;
close all;


L = 50; % ������
M = 8; % ��Ԫ��
k = 2; % ��Դ��
SNR = 5; % �����
d = 0.5; %
t = 0:L-1;
K = sqrt(2*10^(SNR/10)); % �źŷ���
the0 = [40 45];
the1 = the0(1) * pi /180;
the2 = the0(2) * pi /180;

w = [pi/8,pi/6]; % frequency
fai1 = unidrnd(360,1,L)*pi/180;
fai2 = unidrnd(360,1,L)*pi/180;

s1 = K*exp(1i * (w(1) * t + fai1));
s2 = K*exp(1i * (w(2) * t + fai2));
figure;
subplot(1,2,1);
plot(real(K*exp(1i * (w(1) * t))));
title('w = pi/8 �������ź�');
xlabel('������');
ylabel('��ֵ');

subplot(1,2,2);
scatter(t,fai1,'+');
title('�����λ');
xlabel('������');
ylabel('��λ��pi)')

figure;
plot(real(s1));
title('��Դһ�Ĳ���ͼ');
xlabel('������');
ylabel('��ֵ');

Nn = wgn(M,L,0,'complex');
figure;
plot(t,real(Nn(1,:)));
title('��˹������');
xlabel('������');
ylabel('��ֵ');


