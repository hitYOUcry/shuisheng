  %---�ǰ����������Ƶ�ʣ���λ�źŵļ��(2_3)
  %---��ɫ��������Ԥ�׻����ټ��
 
  
 clear all;
 clc
%-----�������ķǰ������Ĺ�������״ȷ��AR�˲�����ϵ����
%-----ȡ������ɢ��ֵ��ʹ֮���ǹ����׵���Ҫ����
  
  fc=37500;                                                  %---�źŵ�Ƶ��
  fs=1000000;                                                 %---����Ƶ��
  Time=0.0001 ;                                            
  N=floor(fs*Time);
  k=0.1;	
  fm=37400;	 
  f0=37600; 
  f=0:fs/2-1;
  G_f=( (fm+k*(f+f0))./(fm^2+(f+f0).^2) + (fm-k*(f-f0))./(fm^2+(f-f0).^2))./pi;  %--����������ɢ��
  G_f=[G_f,fliplr(G_f)];
  G_n=real(ifft(G_f));
  [AR,E]=levinson(G_n,4);
  AR
  
  figure(1)
  plot(G_f,'b'),xlabel('Ƶ��'),ylabel('G_f');
  %title('The shape of Gf')           %����G_f
  hold on
  b=sqrt(E);
  noisefreq=freqz(b,AR,fs,'whole');
  noisePsd=noisefreq.*conj(noisefreq);%�ǰ������Ĺ������ܶ�
  plot(noisePsd,'r:')
 % title('����������ģ�������Ĺ�����') ;
  legend('��������','��ģ��õ�������');
  grid on
  hold off
 
  
  
%-----�ڹ۲�ʱ����ͬ����£��Ƚϲ�ͬ�龯���ʺ�����ȶԼ���Ӱ��-----
 
  alpha=0.1;                                                  %---Ҫ����龯���ʵ�����
 
  SNR=-40:0.1:0;
  length=size(SNR,2); 
  pd=zeros(1,length);
  pa=zeros(1,length);
  for i=1:length
     [pd(:,i),pa(:,i)]=test_signal_3(fc,fs,SNR(i),Time,AR,alpha);              %---Ԥ�׻��Ӻ���
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
  [pd(:,i),pa(:,i)]=test_signal_3(fc,fs,SNR(i),Time,AR,alpha);              %---Ԥ�׻��Ӻ���
  end
  figure(2);
  plot(SNR,pd,'b--');

  
  alpha=0.1;                                                             %---Ҫ����龯���ʵ�����
  pd=zeros(1,length);
  pa=zeros(1,length);
  for i=1:length
  [pd(:,i),pa(:,i)]=test_signal_2(fc,fs,SNR(i),Time,AR,alpha);              %---��Ԥ�׻��Ӻ���
  end
  figure(2);
  plot(SNR,pd,'g')
  hold on;
   figure(3);
  plot(SNR,pa,'g');
  hold on;
  
  alpha=0.01;                                                            %---�龯����alpha��0.001
  pd=zeros(1,length);
  pa=zeros(1,length);
  for i=1:length
  [pd(:,i),pa(:,i)]=test_signal_2(fc,fs,SNR(i),Time,AR,alpha);              %---��Ԥ�׻��Ӻ���
  end
  figure(2);
  plot(SNR,pd,'g--')
  


  
  figure(2)
  grid on;
  title({['�ź�Ƶ����f=',num2str(fc),'�ҹ۲�ʱ��ΪTime=',...
         num2str(Time),'�¹�����������']});
  xlabel('�����SNR');
  ylabel('��ȷ������');
  legend('Ԥ��_�龯����Ϊ0.1','Ԥ��_�龯����Ϊ0.05','��Ԥ��_�龯����Ϊ0.1','��Ԥ��_�龯����Ϊ0.01',2);
  
   figure(3)
  grid on;
  title({['�ź�Ƶ����f=',num2str(fc),'�ҹ۲�ʱ��ΪTime=',...
         num2str(Time),'�¹�����������']});
  xlabel('�����SNR');
  ylabel('�龯����');
  %axis([-20 0   -0.001  0.1]);
  legend('Ԥ��_�龯����Ϊ0.1','Ԥ��_�龯����Ϊ0.05','��Ԥ��_�龯����Ϊ0.1','��Ԥ��_�龯����Ϊ0.01');
  
  
  
  


  
  
  
  
  
