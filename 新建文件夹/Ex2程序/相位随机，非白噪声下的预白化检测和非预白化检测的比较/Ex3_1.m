  %---非白噪声中随机频率，相位信号的检测(2_3)
  %---有色噪声进行预白化，再检测
 
  
 clear all;
 clc
%-----由所给的非白噪声的功率谱形状确定AR滤波器的系数，
%-----取若干离散点值，使之覆盖功率谱的主要部分
  
  fc=37500;                                                  %---信号的频率
  fs=1000000;                                                 %---采样频率
  Time=0.0001 ;                                            
  N=floor(fs*Time);
  k=0.1;	
  fm=37400;	 
  f0=37600; 
  f=0:fs/2-1;
  G_f=( (fm+k*(f+f0))./(fm^2+(f+f0).^2) + (fm-k*(f-f0))./(fm^2+(f-f0).^2))./pi;  %--将功率谱离散化
  G_f=[G_f,fliplr(G_f)];
  G_n=real(ifft(G_f));
  [AR,E]=levinson(G_n,4);
  AR
  
  figure(1)
  plot(G_f,'b'),xlabel('频率'),ylabel('G_f');
  %title('The shape of Gf')           %画出G_f
  hold on
  b=sqrt(E);
  noisefreq=freqz(b,AR,fs,'whole');
  noisePsd=noisefreq.*conj(noisefreq);%非白噪声的功率谱密度
  plot(noisePsd,'r:')
 % title('理论噪声和模拟噪声的功率谱') ;
  legend('理论噪声','建模后得到的噪声');
  grid on
  hold off
 
  
  
%-----在观测时间相同情况下，比较不同虚警概率和信噪比对检测的影响-----
 
  alpha=0.1;                                                  %---要求的虚警概率的上限
 
  SNR=-40:0.1:0;
  length=size(SNR,2); 
  pd=zeros(1,length);
  pa=zeros(1,length);
  for i=1:length
     [pd(:,i),pa(:,i)]=test_signal_3(fc,fs,SNR(i),Time,AR,alpha);              %---预白化子函数
  end
  figure(2);
  plot(SNR,pd,'b');
  hold on;
   figure(3);
  plot(SNR,pa,'r');
  hold on;
 
  alpha=0.01;                                                
  pd=zeros(1,length);
  pa=zeros(1,length);
  for i=1:length
  [pd(:,i),pa(:,i)]=test_signal_3(fc,fs,SNR(i),Time,AR,alpha);              %---预白化子函数
  end
  figure(2);
  plot(SNR,pd,'b--');

  
  alpha=0.1;                                                             %---要求的虚警概率的上限
  pd=zeros(1,length);
  pa=zeros(1,length);
  for i=1:length
  [pd(:,i),pa(:,i)]=test_signal_2(fc,fs,SNR(i),Time,AR,alpha);              %---非预白化子函数
  end
  figure(2);
  plot(SNR,pd,'g')
  hold on;
   figure(3);
  plot(SNR,pa,'g');
  hold on;
  
  alpha=0.01;                                                            %---虚警概率alpha是0.001
  pd=zeros(1,length);
  pa=zeros(1,length);
  for i=1:length
  [pd(:,i),pa(:,i)]=test_signal_2(fc,fs,SNR(i),Time,AR,alpha);              %---非预白化子函数
  end
  figure(2);
  plot(SNR,pd,'g--')
  


  
  figure(2)
  grid on;
  title({['信号频率在f=',num2str(fc),'且观测时间为Time=',...
         num2str(Time),'下工作特性曲线']});
  xlabel('信噪比SNR');
  ylabel('正确检测概率');
  legend('预白_虚警概率为0.1','预白_虚警概率为0.05','非预白_虚警概率为0.1','非预白_虚警概率为0.01',2);
  
   figure(3)
  grid on;
  title({['信号频率在f=',num2str(fc),'且观测时间为Time=',...
         num2str(Time),'下工作特性曲线']});
  xlabel('信噪比SNR');
  ylabel('虚警概率');
  %axis([-20 0   -0.001  0.1]);
  legend('预白_虚警概率为0.1','预白_虚警概率为0.05','非预白_虚警概率为0.1','非预白_虚警概率为0.01');
  
  
  
  


  
  
  
  
  
