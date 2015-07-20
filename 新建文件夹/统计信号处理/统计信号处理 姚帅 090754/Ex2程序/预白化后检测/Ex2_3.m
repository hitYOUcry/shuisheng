  %---�ǰ����������Ƶ�ʣ���λ�źŵļ��(2_3)
  %---��ɫ��������Ԥ�׻����ټ��
  %---�����ۣ�2008-5
  
 clear all;

%-----�������ķǰ������Ĺ�������״ȷ��AR�˲�����ϵ����
%-----ȡ������ɢ��ֵ��ʹ֮���ǹ����׵���Ҫ����
  k=0.1;	
  fm=200;	 
  f0=260; 
  f=0:9999;
  G_f=( (fm+k*(f+f0))./(fm^2+(f+f0).^2) + (fm-k*(f-f0))./(fm^2+(f-f0).^2))./pi;  %--����������ɢ��
  G_n=abs(ifft(G_f,10000));
 % G_n1=G_n(1:5);
  [AR,E]=levinson(G_n,4);
  AR                                                    %---�õ�AR�˲�����ϵ��
 
  fc=240;                                                  %---�źŵ�Ƶ��
  ph=pi/4;
  fs=4000;                                                 %---����Ƶ��
  Time=0.01                                               %---����ʱ��Ϊ��ͬ�����£���ͬ����ȼ�����  
  
  
 %Ԥ�׻�
  a=0.01;                                                 
  pd=zeros(1,20);
  pa=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_3(fc,ph,fs,SNR(i),Time,AR,a);              %---���ù��ܺ�����ʵ��
  end
  figure(1);
  plot(SNR,pd,'b');
  hold on;
  figure(2);
  plot(SNR,pa,'b');
  hold on;
  
  %��Ԥ�׻�
  a=0.01;                                               
  pd=zeros(1,20);
  SNR=-15:4;
  for i=1:20
  pd(:,i)=test_signal_2(fc,ph,fs,SNR(i),Time,AR,a);              %---���ù��ܺ�����ʵ��
  end
  figure(1);
  plot(SNR,pd,'m:')
  figure(3);
  plot(SNR,pa,'m:');
  
 figure(1)
 grid on;
 % title({['�ź�Ƶ����f=',num2str(fc),'���ҹ۲�ʱ��ΪTime=',...
  %       num2str(Time),'�¹�����������']});
  xlabel('�����SNR');
  ylabel('��ȷ������');
  legend('Ԥ��','��Ԥ��');
  
  figure(2)
  grid on;
 % title({['�ź�Ƶ����f=',num2str(fc),'���ҹ۲�ʱ��ΪTime=',...
 %        num2str(Time),'�¹�����������']});
  xlabel('�����SNR');
  ylabel('�龯����');
  %axis([-20 0   -0.001  0.1]);
  legend('Ԥ��','��Ԥ��');
  
  
  
  


  
  
  
  
  
