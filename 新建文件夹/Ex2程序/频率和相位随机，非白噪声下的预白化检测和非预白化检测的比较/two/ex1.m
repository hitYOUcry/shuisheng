 %---在非白噪声情况下，随机频率，随机相位信号，进行5000次实验，进行检测,没有进行预白化的情况
  
  clear all;
  clc

%-----由所给的非白噪声的功率谱形状确定AR滤波器的系数，
%-----取若干离散点值，使之覆盖功率谱的主要部分
  k=0.05;	
  fm=200;	 
  f0=320; 
  fs=4000;
  f=0:fs/2-1;
  G_f=( (fm+k*(f+f0))./(fm^2+(f+f0).^2) + (fm-k*(f-f0))./(fm^2+(f-f0).^2))./pi;  %--将功率谱离散化
  G_f=[G_f,fliplr(G_f)];
  G_n=real(ifft(G_f));
 % G_n1=G_n(1:5);
  [AR,E]=levinson(G_n,4);
  AR                                                    %---得到AR滤波器的系数           
        
    
                                              
 
  fs=4000;                                                 %---采样频率
  SNR=-10:0.5:-0.5;
  Time=0.05;
%-----在观测时间相同情况下，比较不同虚警概率和信噪比对检测的影响-----
  alpha=0.01;
 for i=1:20
 [ pd1(:,i), pa1(:,i)]=n_yubai(fs,SNR(i),Time,AR,alpha);              
 [ pd2(:,i), pa2(:,i)]=yubai(fs,SNR(i),Time,AR,alpha); 
 end
   figure(1)
  plot(SNR,pd1,'r-.')
  hold on
  
  plot(SNR,pd2,'g-.')
  figure(2)
  plot(SNR,pa1,'r-.')
  hold on
  plot(SNR,pa2,'g-.')
  
   alpha=0.02;                                         
  for i=1:20
    [ pd3(:,i), pa3(:,i)]=n_yubai(fs,SNR(i),Time,AR,alpha);             
    [ pd4(:,i), pa4(:,i)]=yubai(fs,SNR(i),Time,AR,alpha); 
  end
  figure(1)
  plot(SNR,pd3,'r')
   plot(SNR,pd4,'g')
  figure(2)
   plot(SNR,pa3,'r')
    plot(SNR,pa4,'g')
 
  
  figure(1)
  grid on
  title('随机频率和随机相位信号在有色噪声下,不同门限下的检测概率');
  xlabel('信噪比SNR');
  ylabel('正确检测概Pd');
  legend('非预白化alpoha=0.01','预白化alplha=0.01','非预白化alpha=0.02','预白化alpha=0.02');
  
  figure(2)
   grid on;
   title('随机频率和随机相位信号在有色噪声下,不同门限下的虚警概率');
  xlabel('信噪比SNR');
  ylabel('虚警概率alpha')
  axis([-11  0  0  0.12]);
   legend('非预白化alpha=0.01','预白化alpha=0.01','非预白化alpha=0.02','预白化alpha=0.02');
  
  
  figure(3)
  plot(SNR,pd1./pa1,'rp')
  hold on
   plot(SNR,pd2./pa2,'g*')
  
  
  
