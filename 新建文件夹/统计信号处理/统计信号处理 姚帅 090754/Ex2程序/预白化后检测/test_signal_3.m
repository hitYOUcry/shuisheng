%---非白噪声中随机频率，相位信号的检测(子函数)

function  Pd=test_signal_3(fc,ph,fs,SNR,Time,AR,a)
   N=floor(fs*Time);                                %---观测时间里的采样点数
   f=fc/fs;                                         %---归一化信号的频率
   A=1;                                             %---固定信号的幅值，信噪比的改变可通过改变噪声的功率来实现
   S=A.*sin(2*pi*f*[1:N]+ph);
   Es = S*S';                                       %---信号能量
   Ps=Es/N;
   Dn = Ps/power(10,SNR/10);                        %---噪声功率
   Z0= 0.5*sqrt(-2*Dn*Es*log(a));                   %---判决门限，由虚警概率上限求出，可在Z0前乘一系数来改变Z0
   
%---------———————用5000次实验来计算检测情况---------——————————
   Pd=0;                                            %---Pd用于存放正确检测概率
   Pa=0;                                            %---Pa用于存放虚警概率
 for k=1:5000
   G_W_noise=wgn(1,N,1,'linear');                  %---产生单位高斯白噪声          
   noise=filter(1,AR,G_W_noise);                   %---产生所要求的谱形状的非白噪声
  
   sum1=0;                                          %---本小段程序是把修改噪声的功率，使之满足SNR
   for i=1:N
       sum1=sum1+noise(i).^2;
   end
   sum1=sum1/N;
   para=Dn./sum1;
   noise=noise*sqrt(para);
  
   x=A.*sin(2*pi*f*[1:N]+ph)+noise;                 %---接收到的有信号的随机过程
   x=filter(AR,1,x);                                %---预白化
   Ss=A.*sin(2*pi*f*[1:N]);                  
   Sc=A.*cos(2*pi*f*[1:N]);
   Gss(k)=x*Ss'; 
   Gcs(k)=x*Sc';
   Zs(k)=sqrt(Gss(k).^2+Gcs(k).^2);               
   if Zs(k)>Z0
       Pd=Pd+1;
   end
   
   noise=filter(AR,1,noise);
   Gsn(k)=noise*Ss';                                 %---接收到的无信号的随机过程
   Gsc(k)=noise*Sc';
   Zn(k)=sqrt(Gsn(k).^2+Gsc(k).^2);
   if Zn(k)>Z0
       Pa=Pa+1;
   end
 end
  Pd=Pd/5000;
  Pa=Pa/5000;
