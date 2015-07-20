%基本ESPRIT算法,第二种方法
clear all;close all;clc;

c=3*10^8;
f=3*10^9;
lamda=c/f;
d=lamda/2;
n=10;
signal_number=3;

thita1=-25;
thita2=30;
thita3=65;

f1=40;
f2=20;
f3=70;

snapshot=1:2000;

S1=4.4*exp(j*2*pi*f1*snapshot/length(snapshot)); 
S2=2.5*exp(j*2*pi*f2*snapshot/length(snapshot)); 
S3=3.6*exp(j*2*pi*f3*snapshot/length(snapshot)); 
S=[S1;S2;S3];

A1=exp(-j*2*pi*d*[0:n-1]*sin(thita1*pi/180)/lamda).';
A2=exp(-j*2*pi*d*[0:n-1]*sin(thita2*pi/180)/lamda).';
A3=exp(-j*2*pi*d*[0:n-1]*sin(thita3*pi/180)/lamda).';
A=[A1,A2,A3];

N=sqrt(1/2)*(randn(n,length(snapshot))+j*randn(n,length(snapshot)));

X=A*S+N;

Rxx=X*X'/length(snapshot);

[V,D]=eig(Rxx);
[Y,I]=sort(diag(D));
Us=V(:,n-signal_number+1:n);
U1=Us(1:n-1,:);
U2=Us(2:n,:);
[p,q]=eig(inv(U1'*U1)*U1'*U2);          %张贤达《矩阵分析与应用》 第528页

for i=1:signal_number;
    alpha(i)=real(asin(-j*(log(q(i,i)))*lamda/(-2*pi*d))*180/pi);
end;
figure(2)
stem(alpha,ones(1,signal_number),'r--');grid;
axis([-90 90 0 2]);
text(alpha(1)-4,1.1,num2str(alpha(1)));
text(alpha(2)-4,1.1,num2str(alpha(2)));
text(alpha(3)-4,1.1,num2str(alpha(3)));
title('ESPRIT算法DOA估计');




