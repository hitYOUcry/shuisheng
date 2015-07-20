 %---在非白噪声情况下，随机频率，随机相位信号，进行5000次实验，进行检测,没有进行预白化的情况
  %---张艳芹，2008-5
  
  clear all;
  clc
   fs=1000000; 
%-----由所给的非白噪声的功率谱形状确定AR滤波器的系数，
%-----取若干离散点值，使之覆盖功率谱的主要部分
k=0.1;	
fm=37400;	 
f0=37600; 
nfft=2^(floor(log2(fs))+1);           %对功率谱进行傅里叶反变换的点数
f=linspace(-fs/2,fs/2,nfft+1);        %f=linspace(0,fs,nfft+1)-fs/2;
f=f(1:end-1);
G_f=( (fm+k*(f+f0))./(fm^2+(f+f0).^2) + (fm-k*(f-f0))./(fm^2+(f-f0).^2))./pi;  %--将功率谱离散化
G_n=[G_f(nfft/2+1:end),G_f(1:nfft/2)];     
% G_n1=G_n(1:5);
rx=ifft(G_n);                          %对功率谱进行傅里叶反变换的自相关序列rx(k)
rxk=rx(1:5);
[AR,E]=levinson(rxk,4);                                                      %---得到AR滤波器的系数
  fc=37500+rand(1)*100;                                                  %---信号的频率               
  %---采样频率
   SNR=-40:0.1:0;
   length=size(SNR,2);
  alpha=0.01;
%-----在观测时间相同情况下，比较不同虚警概率和信噪比对检测的影响-----

 Time=0.0001;
  for i=1:length
        [ pd1(:,i), pa1(:,i)]=test_signal_2(fc,fs,SNR(i),Time,AR,alpha);              %---调用功能函数来实现
  end
%  
%   fs=20000000;
%   for i=1:length
%   [ pd2(:,i), pa2(:,i)]=test_signal_2(fc,fs,SNR(i),Time,AR,alpha);              %---调用功能函数来实现
%   end
%   fs=50000000;
%    for i=1:length
%   [ pd3(:,i), pa3(:,i)]=test_signal_2(fc,fs,SNR(i),Time,AR,alpha);              %---调用功能函数来实现
%    end
 
      figure(1)
    plot(SNR,pd1,'r')
    hold on
    grid on
 %     plot(SNR,pd2,'g')
 %     hold on
 %     plot(SNR,pd3,'b')
  title('有色噪声干扰下不同采样频率下随机相位、频率信号实际虚警概率');
  xlabel('SNR(db)');
  ylabel('系统正确检测概Pd');
%    legend('采样频率为1MHZ','采样频率20MHZ','采样频率50MHZ');
  
     figure(2)
    plot(SNR,pa1,'r')
 %     hold on
 %     grid on
%      plot(SNR,pa2,'g')
%      hold on
 %     plot(SNR,pa3,'b')
  title('有色噪声干扰下不同采样频率下随机相位、频率信号虚警检测概率');
  xlabel('SNR(db)');
  ylabel('系统虚警概率Pf');
 %   legend('采样频率为1MHZ','采样频率20MHZ','采样频率50MHZ');
  
  
