%---非白噪声中随机频率，相位信号的检测(子函数)
%---有色噪声下的检测

function  [Pd,Pa]=test_signal_2(fc,fs,SNR,Time,AR,alpha)
   N=floor(fs*Time);                                %---观测时间里的采样点数
   A=1;                                             %---固定信号的幅值，信噪比的改变可通过改变噪声的功率来实现
   f=fc/fs;
   ph=2*pi*rand;
   S=A.*sin(2*pi*(240/fs)*[1:N]+ph);
   Es = S*S';                                       %---信号能量
   Ps=Es/N;
   Dn = Ps/power(10,SNR/10);                        %---噪声功率
   Z0=1.85* sqrt(-Dn*Es*log(alpha));                   %---判决门限，由虚警概率上限求出，可在Z0前乘一系数来改变Z0
   
%---------———————用5000次实验来计算检测情况---------——————————
   Pd=0;                                            %---Pd用于存放正确检测概率
   Pa=0;                                            %---Pa用于存放虚警概率
 for k=1:500
   G_W_noise=wgn(1,N,1,'linear');                 %---产生单位高斯白噪声          
   noise=filter(1,AR,G_W_noise);                  %---产生所要求的谱形状的非白噪声
  
   sum1=0;                                          %---本小段程序是把修改噪声的功率，使之满足SNR
   for i=1:N
       sum1=sum1+noise(i).^2;
   end
   sum1=sum1/N;
   para=Dn./sum1;
   noise=noise*sqrt(para);
   ph=2*pi*rand;
   x=A.*sin(2*pi*f*[1:N]+ph)+noise;                 %---接收到的有信号的随机过程
   for w=0:10
    
        Ss=A*sin(2*pi*((220+w*4)/fs)*[1:N]);                  
        Sc=A*cos(2*pi*((220+w*4)/fs)*[1:N]);
        Gss(k,w+1)=x*Ss'; 
        Gcs(k,w+1)=x*Sc';
        Z_s(k,w+1)=sqrt(Gss(k,w+1).^2+Gcs(k,w+1).^2);                  %---信号的Z变换
        
       Gsn(k,w+1)=noise*Ss';                                 %---接收到的无信号的随机过程
       Gsc(k,w+1)=noise*Sc';
       Z_n(k,w+1)=sqrt(Gsn(k,w+1).^2+Gsc(k,w+1).^2);
   end
   Zs(k)=max(Z_s(k,:));
   Zn(k)=max(Z_n(k,:));
   if Zs(k)>Z0
       Pd=Pd+1;
   end

   if Zn(k)>Z0
       Pa=Pa+1;
   end
   
 end
  Pd=Pd/500;
  Pa=Pa/500;
 
  
   


