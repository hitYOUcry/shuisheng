%---非白噪声中随机频率，随机相位信号的检测(子函数)
%---有色噪声下的检测

function  [Pd,Pa]=yubai(fs,SNR,Time,AR,alpha)
    N=floor(fs*Time);                                %---观测时间里的采样点数
    A=1;                                             %---固定信号的幅值，信噪比的改变可通过改变噪声的功率来实现
  
   
%---------———————用5000次实验来计算检测情况---------——————————
   Pd=0;                                            %---Pd用于存放正确检测概率
   Pa=0;                                            %---Pa用于存放虚警概率
    w=1:41;
    ph=2*pi*rand;
    S=A.*sin((2*pi*(220+w-1)/fs)'*[1:N]+ph); %每行都是一个时间序列，有21行
   tem=S.*S; 
   for w=1:41
         Es(w,1)=sum(tem(w,:));
    end                                 %---信号能量
    Ps=Es/N;
    Dn= Ps/power(10,SNR/10);                        %---噪声功率，是一个列向量
    Z0=1.1*sqrt(-Dn.*Es*log(alpha));               %---判决门限，由虚警概率上限求出，可在Z0前乘一系数来改变Z0。是一个列向量
   jj=0;
   for k=1:500
      f=220+floor(unifrnd(0,40));
       ph=2*pi*rand;
      s1=A.*sin(2*pi*f/fs*[1:N]+ph);
      Es1=s1*s1';
      Ps1=Es1/N;
      Dn1= Ps1/power(10,SNR/10);  
      G_W_noise=wgn(1,N,1,'linear');                 %---产生单位高斯白噪声    
    
       noise=filter(1,AR,G_W_noise);                  %---产生所要求的谱形状的非白噪声
                                              %---本小段程序是把修改噪声的功率，使之满足SNR
        sum1=noise*noise';
       sum1=sum1/N;
       para=Dn1./sum1;
       noise=noise.*sqrt(para);
       
      
      if (rand>0.5)
             x=s1+noise;  
             ii=1;
             jj=jj+1;
      else
                 
             x=noise;                                        %接收到的只是噪声
            ii=0;
      end
      x=filter(AR,1,x) ;                               %---预白化
      
        w=1:41;
        Ss=A*sin((2*pi*(220+w-1)/fs)'*[1:N]); 
        Sc=A*cos((2*pi*(220+w-1)/fs)'*[1:N]);
         for n=1:41
            Ss(n)=filter(AR,1,Ss(n));
            Sc(n)=filter(AR,1,Sc(n));
         end
        Gss=x*Ss'; 
        Gcs=x*Sc';
         
        for w=1:41
            Z_s(w)=sqrt(Gss(w).^2+Gcs(w).^2);                
            Z_s(w)= Z_s(w)/Z0(w);               %判决量对相应的门限值归一化   
        end
        Zs =max(Z_s);
        if ((Zs>1)&&(ii==1))
              Pd=Pd+1;
        end
       if ((Zs>1)&&(ii==0))
         Pa=Pa+1;
        end
       
    
 end
  Pd=Pd/jj;
  Pa=Pa/(500-jj);
 
  
     

