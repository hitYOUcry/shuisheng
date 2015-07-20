function [color]=color_produce(k,fm,f0,N)
%%%color_produce: 由给定的功率谱G(f)产生色噪声（信号）
%%%色噪声（信号） ：AR 4阶模型
%%%input:          k,fm,f0------决定谱形状（eg k=0.1,fm=200,f0=260-----噪声；k=） 
%%%                N------------输出的噪声（信号）的个数
%%%output:       :产生给定功率谱的色噪声（信号）
% 方法：1.对G(f)进行傅里叶反变换计算得有色噪声的自相关rx,
%       2.对有色噪声color_n进行AR建模：调用函数levinson获得AR模型参数
%       3.高斯白噪声通过滤波器产生需要的有色噪声color_n
%%%by tanjunhong
%%% 6.7,2011, in lab

%%%% test this function   %%%%%%
% % k=0.1;fm=200;f0=260;N=40;       %色噪声功率谱G(f)参数选择
% % k=1;fm=400;f0=360;N=40;         %色信号功率谱G(f)参数选择

%参数
fs=4000;                             %采样频率
p=4;                                 %AR建模的阶数

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  对G(f)进行傅里叶反变换计算得自相关rx  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfft=1024*8;                          %对功率谱进行傅里叶反变换的点数
f=linspace(-fs/2,fs/2,nfft+1);        %f=linspace(0,fs,nfft+1)-fs/2;
f=f(1:end-1);
G=1/pi.*((fm+k*(f+f0))./(fm^2+(f+f0).^2)+(fm-k*(f-f0))./(fm^2+(f-f0).^2))*fs;   % 给定的功率谱G(f)*fs得原模拟信号：[-pi,pi]
G1=[G(nfft/2+1:end),G(1:nfft/2)];                                               %[0,2*pi]

rx=ifft(G1);                         %对功率谱进行傅里叶反变换的自相关序列rx(k)
% disp('autorelation sequence:')
rxk=rx(1:p+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 要求噪声的协方差矩阵C（真实值，不是估计值，估计值由接收的序列进行自相关估计，再得估计的C）
%%%%% 1-----
% rxN=rx(1:N);
% save 'rnk_N'  rxN                  %保存色噪声的真实的自相关序列（N个值）
% C=toeplitz(rxN);
%%%%% 2-----已知前P+1个和AR系数，进行自相关外推
% rxk=[1.8874    1.4124    0.8033    0.3397    0.0196];
% p=4;
% rxN(1:5)=rxk;  a= [-0.9699    0.3081   -0.0510    0.0713]; %滤波器系数
% for k=6:N
%     rxN(k)=-rxN(k-1:-1:k-p)*a';
% end
% C=toeplitz(rxN');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   进行AR建模：调用函数levinson获得AR模型参数   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AR,E]=levinson(rxk,p);     
% disp('estimation of filter parameter')
% AR,E,                  
%%%result
%%%noise  k=0.1;fm=200;f0=260
%%%       rx:1.8874    1.4124    0.8033    0.3397    0.0196
%%%       AR:[1.0000   -0.9699    0.3081   -0.0510    0.0713];E:0.7490
%%%signal k=1;fm=400；f0=360;
%%%       rx:1.8874    1.4124    0.8033    0.3397    0.0196
%%%       AR:[1.0000   -0.9699    0.3081   -0.0510    0.0713];E:0.7490


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   产生需要的色噪声（信号） %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wn=randn(N,1);
color=filter(sqrt(E),AR,wn);



% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%  验证AR建模的正确性   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % %%estimated power destiny
% % w=2*pi/fs.*f;
% % G_est=E./abs(1+AR(2).*exp(-1i*w)+AR(3).*exp(-2*1i*w)+AR(4).*exp(-3*1i*w)+AR(5).*exp(-4*1i*w)).^2;
% % 
% % % % % figure
% % % % % plot(2*f/fs,10*log10(G),'b');
% % % % % grid on;hold on
% % % % % xlabel('normalized frequency  / pi');ylabel('amplitude / dB');
% % % % % title('the true power destiny');
% % % % % plot(2*f/fs,10*log10(GW),'r:');
% % % % % legend('true power destiny','estimated power destiny')
% % figure(2)
% % plot(f,10*log10(G),'b');hold on
% % plot(f,10*log10(G_est),'r:');grid on;
% % legend('true power destiny','estimated power destiny');
% % xlabel(' frequency  / Hz');ylabel('amplitude / dB');
% % title('power destiny : k=0.1 fm=200 f0=260 ');
% % % title('power destiny : k=1 fm=400 f0=360 ');









%附注：计算的有色噪声的自相关函数rx
% % % % %%%----------- method 1---------------------------%%%%%%%%%%%%%%
%%%%%%%%%%%%%-----------matlab符号运算-------------------------- %%%%%%%%%%%%%%
% % w0=2*pi*f0; %模拟角频率
% % wm=2*pi*fm;
% % syms w
% % G=1/pi*((fm+k*(w/2/pi+f0))/(fm^2+(w/2/pi+f0)^2)+(fm-k*(w/2/pi-f0))/(fm^2+(w/2/pi-f0)^2));
% % ff=ifourier(G);
% % syms x
% % ff=simplify(ff);
% % x=0:1/fs:4/fs;
% % ff=exp(pi*x*(520*i + 400)).*(i/10 - heaviside(x).*(i/10 + 1) + 1) + (heaviside(x)*(i/10 - 1) - i/10 + 1)./exp(pi*x*(520*i - 400)) + (2*cos(520*pi*x).*heaviside(x))./exp(400*pi*x) + (sin(520*pi*x).*heaviside(x))./(5*exp(400*pi*x));
% % ff    % 2.0000    1.3987    0.8082    0.3373    0.0211
% % % % %%% %%%-----------method 2-------------------------------------%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%自己计算的自相关序列，对功率谱进行傅里叶反变换得到自相关序列%%%%%%%%%%%%
% rx0=2;
% t=1/fs:1/fs:p/fs;
% rx=2*(cos(w0.*t)+k*sin(w0.*t)).*exp(-wm.*t).*heaviside(t)+2*(cos(w0.*t)-k*sin(w0.*t)).*exp(wm.*t).*heaviside(-t); 
% rx=[rx0,rx];% 2.0000    1.3987    0.8082    0.3373    0.0211   

%%%%%%%%%%%%%%using method 1 and 2 ,estimated  parameter of AR filter:
% %  AR=[1.0000   -0.7963    0.0976   0.0507 0.0489 ];  E=0.9832
%%%%%%%%%%%%%%%estimated power destiny is not good enough





