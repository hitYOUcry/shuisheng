clear all;
fs = 4000;                          %采样频率(Hz)
fc = 250;                           %信号频率(Hz)
a = 1/10;                           %虚警概率上限
SNR=[-20:1:5];                             %信噪比
Tlen=10/1000;                       %观测时间T
N_Exp=5000;                         %实验次数

figure;
while(a>=1/1000)
    for(i=1:1:length(SNR))
        Pd(i)=WGN_phasetest(fs,fc,a,SNR(i),Tlen,N_Exp)
        i
    end
    if(a==1/10)
         plot(SNR,Pd,'r');
         title({'接收机的工作特性曲线',...
              ['采样频率fs(Hz) = ',num2str(fs),', 信号频率fc(Hz) = ',num2str(fc),...
              ', 观测时间Tlen(ms) = ',num2str(Tlen*1000)]});
          hold on;
      end
      if(a==1/100)   
          plot(SNR,Pd,'g');
          hold on;
      end
      if(a==1/1000)   
          plot(SNR,Pd,'b');
      end
      a=a/10;
end