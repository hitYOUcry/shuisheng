 %--- 要求虚警概率不大于 1% ，在高斯白噪声和信号频率确定相位随机的假设下设计门限Z0，计算正确检测概率

  
  clear all;
 
  fc=240;                                                  %---信号的频率
  fs=4000;                                                 %---采样频率
  alpha=0.01;                                             %---要求的虚警概率的上限 
  
  
%-----在不同观测时间相同情况下，比较检测概率-----
 
  Time=0.01;                                             %---观测时间为0.01s                                                
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,fs,SNR(i),Time,alpha);            
  end
  figure(2);
  plot(SNR,pd,'b');
  hold on;
  
  Time=0.05;                                                    %---观测时间为0.05s                                
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,fs,SNR(i),Time,alpha);             
  end
  plot(SNR,pd,'r-.')
  hold on;
  
  Time=0.1;                                                  %---观测时间为0.1s
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,fs,SNR(i),Time,alpha);             
  end
  plot(SNR,pd,'g:')
  
  grid on;
  title({['信号频率在f=',num2str(fc),'虚警概率alpha=',...
         num2str(alpha),'下工作特性曲线']});
  xlabel('信噪比SNR');
  ylabel('正确检测概率Pd');
  legend('Time=0.01','Time=0.05','Time=0.1');
  
  
  
  