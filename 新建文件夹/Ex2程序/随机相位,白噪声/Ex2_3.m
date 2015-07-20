  clc;
  clear all;
  fc=240;                                                  %---信号的频率
  fs=4000;                                                 %---采样频率
  alpha=0.01;                                              %---要求的虚警概率的上限
  SNR=-15:4;
 
  f=fc/fs;                                        
  A=1;
  ph=2*pi*rand;
  
  figure(3);
  Time=0.01 ; 
  Z0=zeros(1,20);
  N=floor(fs*Time);                                 %---观测时间里的采样点数
     
      S=A.*sin(2*pi*f*[1:N]+ph);
      Es = S*S';                                        %---信号能量
      Ps=Es/N;
      Dn = Ps./(10.^(SNR/10));                         %---噪声功率
      Z0=sqrt(-Dn*Es*log(alpha));                     %---判决门限             
     
      plot(SNR,Z0,'b');
      hold on;
     
      Time=0.05 ; 
      N=floor(fs*Time);                                 %---观测时间里的采样点数
      S=A.*sin(2*pi*f*[1:N]+ph);
      Es = S*S';                                        %---信号能量
      Ps=Es/N;
      Dn = Ps./(10.^(SNR/10));                         %---噪声功率
      Z0=sqrt(-Dn*Es*log(alpha));                     %---判决门限             
      plot(SNR,Z0,'r');
     
    
      Time=0.1;
      N=floor(fs*Time);                                 %---观测时间是0.1s
      S=A.*sin(2*pi*f*[1:N]+ph);
      Es = S*S';                                        %---信号能量
      Ps=Es/N;
      Dn = Ps./(10.^(SNR/10));                         %---噪声功率
      Z0=sqrt(-Dn*Es*log(alpha));                     %---判决门限             

     plot(SNR,Z0,'g');
     
  grid on;
  title({['信号频率在f=',num2str(fc),'虚警概率alpha=',...
         num2str(alpha),'下不同观测时间下的门限值']});
  xlabel('信噪比SNR');
  ylabel('门限');
  legend('Time=0.01','Time=0.05','Time=0.1');
  
   
 