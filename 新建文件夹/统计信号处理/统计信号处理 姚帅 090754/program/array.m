%白噪声时的波束形成
clear
clc
%改进
%只取有数据的一段相加，补零的不要

theta=120;       %目标方位

M=40;           %阵元个数
d=2;            %阵元间距
c=1500;         %声速

distance(1)=10000; %目标与第一个基元的位置
for k=1:M-1
    distance(k+1)=sqrt(distance(1)^2+(k*d)^2-2*distance(1)*k*d*cos(theta*pi/180));    %余弦定理
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

Ps=0.5*A^2;                 %信号功率
SNR=-5;             %信噪比
Pn=Ps/10^(SNR/10);  %噪声功率
sigma2=Pn;
sigma=sqrt(sigma2);


for i=1:M
% tao(i)=-(i-1)*d*cos(theta*pi/180)/c;  %各阵元的延时量
x(i,:)=A*sin(2*pi*f*(n-tao(i)*fs)/fs+fai)+sigma*randn(1,N); %产生阵元数据，用tao做延时,每一路的噪声单独产生
end

win=hanning(M+1)';  %矩形窗效果不好，用hanning窗
k=0:M;

for theta_swap=1:180
    i=1:M;
    tao=-(i-1)*d*cos(theta_swap*pi/180)/c;  %各阵元的延时量    
    %延时滤波器
    tao_n=tao*fs;
    tao_zheng=round(tao_n);     %整数延时
    tao_ling=tao_n-tao_zheng;   %小数延时，tao_ling的值都在-0.5~0.5之间
    y=0;
    for j=2:M
        %第一路可以不要，想要的话可以单独加x(1,:)
        %整数延时
        s_temp(j,:)=x(j,:);        
        if tao_zheng(j)<=-1
             s_after(j,:)=[s_temp(end-abs(tao_zheng(j))+1:end) s_temp(j,1:end-abs(tao_zheng(j)))];    %向后整数延时
        end
        if tao_zheng(j)>=1
            s_after(j,:)=[s_temp(j,1+tao_zheng(j):end) s_temp(1:tao_zheng(j))];    %向前整数延时
        end
        if tao_zheng(j)==0
            s_after(j,:)=s_temp(j,:);    %不做整数延时
        end
        
%         %小数延时
%         hd(j,:)=sin(k*pi-tao_ling(j)*pi)./(k*pi-tao_ling(j)*pi);
%         g(j,:)=hd(j,:).*win;
%         h(j,:)=g(j,:)/sum(g(j,:));
%         s_temp(j,:)=conv(x(j,:),h(j,:));        %小数延时
        y=y+s_after(j,:);   %预形成波束
    end
    P(theta_swap)=y*y';
end
% subplot(2,2,2)
hold on;
plot(P);
title('theta为120度时的波束图');
xlabel('theta');
ylabel('P(theta)');

