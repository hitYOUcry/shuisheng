  
  %--- Ҫ���龯���ʲ����� 1% ���ڸ�˹���������ź�Ƶ��ȷ��,��λ����ļ������������Z0��������ȷ������
 
  
  clear al
  fc=240;                                                  %---�źŵ�Ƶ��
  fs=4000;                                                 %---����Ƶ��
  Time=0.01;                                               %---����ʱ��Ϊ��ͬ�����£���ͬ����ȼ�����  
  
  
%-----�ڹ۲�ʱ����ͬ����£��Ƚ��ڲ�ͬ�龯�����£�����Ⱥ���ȷ�����ʵĹ�ϵ-----
 
  alpha=0.01;                                                  %---�龯���ʵ�������0.01
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,fs,SNR(i),Time,alpha);             
  end
  figure(1);
  plot(SNR,pd,'b');
  hold on;
  
  alpha=0.005;                                                  %---�龯���ʵ�������0.005
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,fs,SNR(i),Time,alpha);             
  end
  plot(SNR,pd,'r-.')
  hold on;
  
  alpha=0.001;                                                  %---�龯���ʵ�������0.001
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,fs,SNR(i),Time,alpha);             
  end
  plot(SNR,pd,'m:')
  
  grid on;
  title({['�ź�Ƶ����f=',num2str(fc),'�ҹ۲�ʱ��ΪTime=',...
         num2str(Time),'�¹�����������']});
  xlabel('�����SNR');
  ylabel('��ȷ������Pd');
  legend('alpah=0.01','alpha=0.005','alpha=0.001');
  
  
  
  
  
  
