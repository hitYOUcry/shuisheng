%---非白噪声中随机频率，相位信号的检测(子函数)

function [Pd,Pa]=test_signal_3(fc,fs,SNR,Time,AR,alpha)
format long;
   N=floor(fs*Time);                                %---观测时间里的采样点数
   f=fc/fs;                                         %---归一化信号的频率
   A=1;             %---固定信号的幅值，信噪比的改变可通过改变噪声的功率来实现
   
   ph=2*pi*rand;
   S=A.*sin(2*pi*f*[1:N]+ph);
   Es = S*S';                                       %---信号能量
   Ps=Es/N;
   Dn = Ps/power(10,SNR/10);                        %---噪声功率
   Z0= 0.24*sqrt(-Dn*Es*log(alpha));                   %---判决门限，由虚警概率上限求出，可在Z0前乘一系数来改变Z0
   
%---------———————用5000次实验来计算检测情况---------——————————
   Pd=0;                                            %---Pd用于存放正确检测概率
   Pa=0;                                            %---Pa用于存放虚警概率
   jj=0;
 for k=1:500
   G_W_noise=wgn(1,N,1,'linear');                  %---产生单位高斯白噪声          
   noise=filter(1,AR,G_W_noise);                   %---产生所要求的谱形状的非白噪声
   ph=2*pi*rand;
   sum1=0;                                          %---本小段程序是把修改噪声的功率，使之满足SNR
   for i=1:N
       sum1=sum1+noise(i).^2;
   end
   sum1=sum1/N;
   para=Dn./sum1;
   noise=noise*sqrt(para);
   ph=2*pi*rand;
   
  
   
   if(rand>0.5)
       x=A.*sin(2*pi*f*[1:N]+ph)+noise;                 %---接收到的有信号的随机过程
       x=filter(AR,1,x);                                %---预白化
       ii=1;
       jj=jj+1;
   else 
        x=noise;                 %---接收到的有信号的随机过程
        x=filter(AR,1,x);  
        ii=0;
   end   
   

   Ss=A.*sin(2*pi*f*[1:N]);   
   Ss=filter(AR,1,Ss);
   Sc=A.*cos(2*pi*f*[1:N]);
   Sc=filter(AR,1,Sc);
   Gss(k)=x*Ss'; 
   Gcs(k)=x*Sc';
   Zs(k)=sqrt(Gss(k).^2+Gcs(k).^2);                  %---信号的Z变换
   if ((Zs(k)>Z0)&&(ii==1))
       Pd=Pd+1;
   end
   if ((Zs(k)>Z0)&&(ii==0))
  
       Pa=Pa+1;
   end
 end
  Pd=Pd/jj;
  Pa=Pa/(500-jj);
 

