%TLS-ESPRIT算法

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

% SNR1=10; 
% SNR2=20; 
% SNR3=30; 

snapshot=1:2000;

S1=4.4*exp(j*2*pi*f1*snapshot/length(snapshot)); 
S2=2.5*exp(j*2*pi*f2*snapshot/length(snapshot)); 
S3=3.6*exp(j*2*pi*f3*snapshot/length(snapshot)); 
S=[S1;S2;S3];

A1=exp(-j*2*pi*d*sin(thita1*pi/180)*[0:n-2]/lamda).';
A2=exp(-j*2*pi*d*sin(thita2*pi/180)*[0:n-2]/lamda).';
A3=exp(-j*2*pi*d*sin(thita3*pi/180)*[0:n-2]/lamda).';
A=[A1,A2,A3];          %子阵1

B1=exp(-j*2*pi*d*sin(thita1*pi/180)*[1:n-1]/lamda).';
B2=exp(-j*2*pi*d*sin(thita2*pi/180)*[1:n-1]/lamda).';
B3=exp(-j*2*pi*d*sin(thita3*pi/180)*[1:n-1]/lamda).';
B=[B1,B2,B3];          %子阵2

N=sqrt(1/2)*(randn(n-1,length(snapshot))+j*randn(n-1,length(snapshot)));

X=A*S+N;
Y=B*S+N;

Rxx=X*X'/length(snapshot);
Rxy=X*Y'/length(snapshot);
[V,D]=eig(Rxx);
s=0;
for i=1:n-1-signal_number;
    s=D(i,i)+s;
end;
delta=s/(n-1-signal_number);

Cxx=Rxx-delta*eye(size(Rxx));
Cxy=Rxy-delta*diag([ones(1,n-1-1)],-1);

[U,S,V]=svd(Cxx);
U1=U(:,1:signal_number);
V1=V(:,1:signal_number);
S1=S(1:signal_number,1:signal_number)
[p,q]=eig(S1,U1'*Cxy*V1);

for i=1:signal_number;
    alpha(i)=real(asin(-j*(log(q(i,i)))*lamda/(-2*pi*d))*180/pi);
end;

stem(alpha,ones(1,signal_number),'r--');grid;
axis([-90 90 0 2]);
text(alpha(1)-4,1.1,num2str(alpha(1)));
text(alpha(2)-4,1.1,num2str(alpha(2)));
text(alpha(3)-4,1.1,num2str(alpha(3)));
title('ESPRIT算法DOA估计');





