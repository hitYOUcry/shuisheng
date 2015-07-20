  %---非白噪声中随机频率，相位信号的检测(2_3)
  %---有色噪声进行预白化，再检测
  %---张艳芹，2008-5
  
 clear all;

%-----由所给的非白噪声的功率谱形状确定AR滤波器的系数，
%-----取若干离散点值，使之覆盖功率谱的主要部分
  k=0.1;	
  fm=200;	 
  f0=260; 
  f=0:9999;
  G_f=( (fm+k*(f+f0))./(fm^2+(f+f0).^2) + (fm-k*(f-f0))./(fm^2+(f-f0).^2))./pi;  %--将功率谱离散化
  G_n=abs(ifft(G_f,10000));
 % G_n1=G_n(1:5);
  [AR,E]=levinson(G_n,4);
  AR                                                    %---得到AR滤波器的系数
 
  fc=240;                                                  %---信号的频率
  ph=pi/4;
  fs=4000;                                                 %---采样频率
  Time=0.01                                               %---计算时间为相同条件下，不同信噪比检测情况  
  
  
 %预白化
  a=0.01;                                                 
  pd=zeros(1,20);
  pa=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_3(fc,ph,fs,SNR(i),Time,AR,a);              %---调用功能函数来实现
  end
  figure(1);
  plot(SNR,pd,'b');
  hold on;
  figure(2);
  plot(SNR,pa,'b');
  hold on;
  
  %非预白化
  a=0.01;                                               
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_2(fc,ph,fs,SNR(i),Time,AR,a);              %---调用功能函数来实现
  end
  figure(1);
  plot(SNR,pd,'m:')
  figure(3);
  plot(SNR,pa,'m:');
  
 figure(1)
 grid on;
 % title({['信号频率在f=',num2str(fc),'秒且观测时间为Time=',...
  %       num2str(Time),'下工作特性曲线']});
  xlabel('信噪比SNR');
  ylabel('正确检测概率');
  legend('预白','非预白');
  
  figure(2)
  grid on;
 % title({['信号频率在f=',num2str(fc),'秒且观测时间为Time=',...
 %        num2str(Time),'下工作特性曲线']});
  xlabel('信噪比SNR');
  ylabel('虚警概率');
  %axis([-20 0   -0.001  0.1]);
  legend('预白','非预白');
  
  
  
  


  
  
  
  
  
