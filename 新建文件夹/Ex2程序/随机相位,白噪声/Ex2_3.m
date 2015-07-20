  clc;
  clear all;
  fc=240;                                                  %---�źŵ�Ƶ��
  fs=4000;                                                 %---����Ƶ��
  alpha=0.01;                                              %---Ҫ����龯���ʵ�����
  SNR=-15:4;
 
  f=fc/fs;                                        
  A=1;
  ph=2*pi*rand;
  
  figure(3);
  Time=0.01 ; 
  Z0=zeros(1,20);
  N=floor(fs*Time);                                 %---�۲�ʱ����Ĳ�������
     
      S=A.*sin(2*pi*f*[1:N]+ph);
      Es = S*S';                                        %---�ź�����
      Ps=Es/N;
      Dn = Ps./(10.^(SNR/10));                         %---��������
      Z0=sqrt(-Dn*Es*log(alpha));                     %---�о�����             
     
      plot(SNR,Z0,'b');
      hold on;
     
      Time=0.05 ; 
      N=floor(fs*Time);                                 %---�۲�ʱ����Ĳ�������
      S=A.*sin(2*pi*f*[1:N]+ph);
      Es = S*S';                                        %---�ź�����
      Ps=Es/N;
      Dn = Ps./(10.^(SNR/10));                         %---��������
      Z0=sqrt(-Dn*Es*log(alpha));                     %---�о�����             
      plot(SNR,Z0,'r');
     
    
      Time=0.1;
      N=floor(fs*Time);                                 %---�۲�ʱ����0.1s
      S=A.*sin(2*pi*f*[1:N]+ph);
      Es = S*S';                                        %---�ź�����
      Ps=Es/N;
      Dn = Ps./(10.^(SNR/10));                         %---��������
      Z0=sqrt(-Dn*Es*log(alpha));                     %---�о�����             

     plot(SNR,Z0,'g');
     
  grid on;
  title({['�ź�Ƶ����f=',num2str(fc),'�龯����alpha=',...
         num2str(alpha),'�²�ͬ�۲�ʱ���µ�����ֵ']});
  xlabel('�����SNR');
  ylabel('����');
  legend('Time=0.01','Time=0.05','Time=0.1');
  
   
 