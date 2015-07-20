clear all;
fs = 4000;                          %����Ƶ��(Hz)
fc = 250;                           %�ź�Ƶ��(Hz)
a = 1/10;                           %�龯��������
SNR=[-20:1:5];                             %�����
Tlen=10/1000;                       %�۲�ʱ��T
N_Exp=5000;                         %ʵ�����

figure;
while(a>=1/1000)
    for(i=1:1:length(SNR))
        Pd(i)=WGN_phasetest(fs,fc,a,SNR(i),Tlen,N_Exp)
        i
    end
    if(a==1/10)
         plot(SNR,Pd,'r');
         title({'���ջ��Ĺ�����������',...
              ['����Ƶ��fs(Hz) = ',num2str(fs),', �ź�Ƶ��fc(Hz) = ',num2str(fc),...
              ', �۲�ʱ��Tlen(ms) = ',num2str(Tlen*1000)]});
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