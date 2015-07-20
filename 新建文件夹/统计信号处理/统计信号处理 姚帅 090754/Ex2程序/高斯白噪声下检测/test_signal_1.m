
%---高斯白噪声下的检测(子函数)

function Pd=test_signal_1(fc,ph,fs,SNR,Time,a)
   N=floor(fs*Time);                                 %---观测时间里的采样点数
   f=fc/fs;                                          %---归一化信号的频率
   A=1;                                           
   S=A.*sin(2*pi*f*[1:N]+ph);
   Es = S*S';                                        %---信号能量
   Ps=Es/N;
   Dn = Ps/power(10,SNR/10);                         %---噪声功率
   Z0=1.1*sqrt(-2*Dn*Es*log(a));                     %---判决门限
 

%---------用5000次实验来计算检测情况---------
   Pd=0;                                             %---Pd用于存放正确检测概率
   Pa=0;                                             %---Pa用于存放虚警概率
 for k=1:5000
   noise=wgn(1,N,Dn,'linear');                       %---产生功率为Dn(watt)的高斯白噪声          
   x=A.*sin(2*pi*f*[1:N]+ph)+noise;                  %---接收到的有信号的随机过程
   Ss=A.*sin(2*pi*f*[1:N]);                          %---确定信号的频率   
   Sc=A.*cos(2*pi*f*[1:N]);
   Gss(k)=x*Ss'; 
   Gcs(k)=x*Sc';
   Zs(k)=sqrt(Gss(k).^2+Gcs(k).^2);                  %---信号的Z变换
   if Zs(k)>Z0
       Pd=Pd+1;
   end

   Gsn(k)=noise*Ss';                                %---接收到无信号的随机过程
   Gsc(k)=noise*Sc';
   Zn(k)=sqrt(Gsn(k).^2+Gsc(k).^2);
   if Zn(k)>Z0
       Pa=Pa+1;
   end
 end
  Pd=Pd/5000;
  Pa=Pa/5000;
  
  
  
   


