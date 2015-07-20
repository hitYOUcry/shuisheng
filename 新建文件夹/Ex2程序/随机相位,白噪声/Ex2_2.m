 %--- Ҫ���龯���ʲ����� 1% ���ڸ�˹���������ź�Ƶ��ȷ����λ����ļ������������Z0��������ȷ������

  
  clear all;
 
  fc=240;                                                  %---�źŵ�Ƶ��
  fs=4000;                                                 %---����Ƶ��
  alpha=0.01;                                             %---Ҫ����龯���ʵ����� 
  
  
%-----�ڲ�ͬ�۲�ʱ����ͬ����£��Ƚϼ�����-----
 
  Time=0.01;                                             %---�۲�ʱ��Ϊ0.01s                                                
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,fs,SNR(i),Time,alpha);            
  end
  figure(2);
  plot(SNR,pd,'b');
  hold on;
  
  Time=0.05;                                                    %---�۲�ʱ��Ϊ0.05s                                
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,fs,SNR(i),Time,alpha);             
  end
  plot(SNR,pd,'r-.')
  hold on;
  
  Time=0.1;                                                  %---�۲�ʱ��Ϊ0.1s
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,fs,SNR(i),Time,alpha);             
  end
  plot(SNR,pd,'g:')
  
  grid on;
  title({['�ź�Ƶ����f=',num2str(fc),'�龯����alpha=',...
         num2str(alpha),'�¹�����������']});
  xlabel('�����SNR');
  ylabel('��ȷ������Pd');
  legend('Time=0.01','Time=0.05','Time=0.1');
  
  
  
  