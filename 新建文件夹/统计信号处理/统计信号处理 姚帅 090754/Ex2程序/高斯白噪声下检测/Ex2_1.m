  
  %--- 要求虚警概率不大于 1% ，在高斯白噪声和信号频率确定的假设下设计门限Z0，计算正确检测概率
  %---张艳芹，2008-5
  
  clear all;
  fc=240;                                                  %---信号的频率
  ph=pi/4;
  fs=4000;                                                 %---采样频率
  Time=0.01;                                               %---计算时间为相同条件下，不同信噪比检测情况  
  
  
%-----在观测时间相同情况下，比较不同虚警概率和信噪比对检测的影响-----
 
  a=0.01;                                                  %---要求的虚警概率的上限
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,ph,fs,SNR(i),Time,a);              %---调用功能函数来实现
  end
  figure(1);
  plot(SNR,pd,'b');
  hold on;
  
  a=0.005;                                                  %---要求的虚警概率的上限
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,ph,fs,SNR(i),Time,a);              %---调用功能函数来实现
  end
  plot(SNR,pd,'r-.')
  hold on;
  
  a=0.001;                                                  %---要求的虚警概率的上限
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,ph,fs,SNR(i),Time,a);              %---调用功能函数来实现
  end
  plot(SNR,pd,'m:')
  
  grid on;
  %title({['信号频率在f=',num2str(fc),'且观测时间为Time=',...
   %      num2str(Time),'下工作特性曲线']});
  xlabel('信噪比SNR');
  ylabel('正确检测概率');
  legend('alpha=0.01','alpha=0.005','alpha=0.001');
  
  
  
  
  
  
