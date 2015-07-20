%������ʱ�Ĳ����γ�
clear
clc
%�Ľ�
%ֻȡ�����ݵ�һ����ӣ�����Ĳ�Ҫ

theta=120;       %Ŀ�귽λ

M=40;           %��Ԫ����
d=2;            %��Ԫ���
c=1500;         %����

distance(1)=10000; %Ŀ�����һ����Ԫ��λ��
for k=1:M-1
    distance(k+1)=sqrt(distance(1)^2+(k*d)^2-2*distance(1)*k*d*cos(theta*pi/180));    %���Ҷ���
end
tao(1)=0;
tao(2:M)=-(distance(1)-distance(2:M))/c;

A=1;
f=200;
fs=4000;
fai=pi/4;
T=1;
N=fs*T;
n=1:N;
s=A*sin(2*pi*f*n/fs+fai);

Ps=0.5*A^2;                 %�źŹ���
SNR=-5;             %�����
Pn=Ps/10^(SNR/10);  %��������
sigma2=Pn;
sigma=sqrt(sigma2);


for i=1:M
% tao(i)=-(i-1)*d*cos(theta*pi/180)/c;  %����Ԫ����ʱ��
x(i,:)=A*sin(2*pi*f*(n-tao(i)*fs)/fs+fai)+sigma*randn(1,N); %������Ԫ���ݣ���tao����ʱ,ÿһ·��������������
end

win=hanning(M+1)';  %���δ�Ч�����ã���hanning��
k=0:M;

for theta_swap=1:180
    i=1:M;
    tao=-(i-1)*d*cos(theta_swap*pi/180)/c;  %����Ԫ����ʱ��    
    %��ʱ�˲���
    tao_n=tao*fs;
    tao_zheng=round(tao_n);     %������ʱ
    tao_ling=tao_n-tao_zheng;   %С����ʱ��tao_ling��ֵ����-0.5~0.5֮��
    y=0;
    for j=2:M
        %��һ·���Բ�Ҫ����Ҫ�Ļ����Ե�����x(1,:)
        %������ʱ
        s_temp(j,:)=x(j,:);        
        if tao_zheng(j)<=-1
             s_after(j,:)=[s_temp(end-abs(tao_zheng(j))+1:end) s_temp(j,1:end-abs(tao_zheng(j)))];    %���������ʱ
        end
        if tao_zheng(j)>=1
            s_after(j,:)=[s_temp(j,1+tao_zheng(j):end) s_temp(1:tao_zheng(j))];    %��ǰ������ʱ
        end
        if tao_zheng(j)==0
            s_after(j,:)=s_temp(j,:);    %����������ʱ
        end
        
%         %С����ʱ
%         hd(j,:)=sin(k*pi-tao_ling(j)*pi)./(k*pi-tao_ling(j)*pi);
%         g(j,:)=hd(j,:).*win;
%         h(j,:)=g(j,:)/sum(g(j,:));
%         s_temp(j,:)=conv(x(j,:),h(j,:));        %С����ʱ
        y=y+s_after(j,:);   %Ԥ�γɲ���
    end
    P(theta_swap)=y*y';
end
% subplot(2,2,2)
hold on;
plot(P);
title('thetaΪ120��ʱ�Ĳ���ͼ');
xlabel('theta');
ylabel('P(theta)');

