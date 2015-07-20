  
  %--- 要求虚警概率不大于 1% ，在高斯白噪声和信号频率确定,相位随机的假设下设计门限Z0，计算正确检测概率
 
  
  clear al
  fc=240;                                                  %---信号的频率
  fs=4000;                                                 %---采样频率
  Time=0.01;                                               %---计算时间为相同条件下，不同信噪比检测情况  
  
  
%-----在观测时间相同情况下，比较在不同虚警概率下，信噪比和正确检测概率的关系-----
 
  alpha=0.01;                                                  %---虚警概率的上限是0.01
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,fs,SNR(i),Time,alpha);             
  end
  figure(1);
  plot(SNR,pd,'b');
  hold on;
  
  alpha=0.005;                                                  %---虚警概率的上限是0.005
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,fs,SNR(i),Time,alpha);             
  end
  plot(SNR,pd,'r-.')
  hold on;
  
  alpha=0.001;                                                  %---虚警概率的上限是0.001
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,fs,SNR(i),Time,alpha);             
  end
  plot(SNR,pd,'m:')
  
  grid on;
  title({['信号频率在f=',num2str(fc),'且观测时间为Time=',...
         num2str(Time),'下工作特性曲线']});
  xlabel('信噪比SNR');
  ylabel('正确检测概率Pd');
  legend('alpah=0.01','alpha=0.005','alpha=0.001');
  
  
  
  
  
  
