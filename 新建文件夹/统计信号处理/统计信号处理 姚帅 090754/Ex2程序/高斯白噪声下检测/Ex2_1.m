  
  %--- Ҫ���龯���ʲ����� 1% ���ڸ�˹���������ź�Ƶ��ȷ���ļ������������Z0��������ȷ������
  %---�����ۣ�2008-5
  
  clear all;
  fc=240;                                                  %---�źŵ�Ƶ��
  ph=pi/4;
  fs=4000;                                                 %---����Ƶ��
  Time=0.01;                                               %---����ʱ��Ϊ��ͬ�����£���ͬ����ȼ�����  
  
  
%-----�ڹ۲�ʱ����ͬ����£��Ƚϲ�ͬ�龯���ʺ�����ȶԼ���Ӱ��-----
 
  a=0.01;                                                  %---Ҫ����龯���ʵ�����
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,ph,fs,SNR(i),Time,a);              %---���ù��ܺ�����ʵ��
  end
  figure(1);
  plot(SNR,pd,'b');
  hold on;
  
  a=0.005;                                                  %---Ҫ����龯���ʵ�����
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,ph,fs,SNR(i),Time,a);              %---���ù��ܺ�����ʵ��
  end
  plot(SNR,pd,'r-.')
  hold on;
  
  a=0.001;                                                  %---Ҫ����龯���ʵ�����
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_1(fc,ph,fs,SNR(i),Time,a);              %---���ù��ܺ�����ʵ��
  end
  plot(SNR,pd,'m:')
  
  grid on;
  %title({['�ź�Ƶ����f=',num2str(fc),'�ҹ۲�ʱ��ΪTime=',...
   %      num2str(Time),'�¹�����������']});
  xlabel('�����SNR');
  ylabel('��ȷ������');
  legend('alpha=0.01','alpha=0.005','alpha=0.001');
  
  
  
  
  
  
