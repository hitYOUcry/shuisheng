 %---�ڷǰ���������£����Ƶ�ʣ������λ�źţ�����5000��ʵ�飬���м��,û�н���Ԥ�׻������
  
  clear all;
  clc

%-----�������ķǰ������Ĺ�������״ȷ��AR�˲�����ϵ����
%-----ȡ������ɢ��ֵ��ʹ֮���ǹ����׵���Ҫ����
  k=0.05;	
  fm=200;	 
  f0=320; 
  fs=4000;
  f=0:fs/2-1;
  G_f=( (fm+k*(f+f0))./(fm^2+(f+f0).^2) + (fm-k*(f-f0))./(fm^2+(f-f0).^2))./pi;  %--����������ɢ��
  G_f=[G_f,fliplr(G_f)];
  G_n=real(ifft(G_f));
 % G_n1=G_n(1:5);
  [AR,E]=levinson(G_n,4);
  AR                                                    %---�õ�AR�˲�����ϵ��           
        
    
                                              
 
  fs=4000;                                                 %---����Ƶ��
  SNR=-10:0.5:-0.5;
  Time=0.05;
%-----�ڹ۲�ʱ����ͬ����£��Ƚϲ�ͬ�龯���ʺ�����ȶԼ���Ӱ��-----
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
  title('���Ƶ�ʺ������λ�ź�����ɫ������,��ͬ�����µļ�����');
  xlabel('�����SNR');
  ylabel('��ȷ����Pd');
  legend('��Ԥ�׻�alpoha=0.01','Ԥ�׻�alplha=0.01','��Ԥ�׻�alpha=0.02','Ԥ�׻�alpha=0.02');
  
  figure(2)
   grid on;
   title('���Ƶ�ʺ������λ�ź�����ɫ������,��ͬ�����µ��龯����');
  xlabel('�����SNR');
  ylabel('�龯����alpha')
  axis([-11  0  0  0.12]);
   legend('��Ԥ�׻�alpha=0.01','Ԥ�׻�alpha=0.01','��Ԥ�׻�alpha=0.02','Ԥ�׻�alpha=0.02');
  
  
  figure(3)
  plot(SNR,pd1./pa1,'rp')
  hold on
   plot(SNR,pd2./pa2,'g*')
  
  
  
