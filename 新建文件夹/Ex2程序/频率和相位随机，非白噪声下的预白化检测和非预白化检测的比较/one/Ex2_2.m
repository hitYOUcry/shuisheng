 %---�ڷǰ���������£����Ƶ�ʣ������λ�źţ�����5000��ʵ�飬���м��,û�н���Ԥ�׻������
  %---�����ۣ�2008-5
  
  clear all;
  clc
   fs=1000000; 
%-----�������ķǰ������Ĺ�������״ȷ��AR�˲�����ϵ����
%-----ȡ������ɢ��ֵ��ʹ֮���ǹ����׵���Ҫ����
k=0.1;	
fm=37400;	 
f0=37600; 
nfft=2^(floor(log2(fs))+1);           %�Թ����׽��и���Ҷ���任�ĵ���
f=linspace(-fs/2,fs/2,nfft+1);        %f=linspace(0,fs,nfft+1)-fs/2;
f=f(1:end-1);
G_f=( (fm+k*(f+f0))./(fm^2+(f+f0).^2) + (fm-k*(f-f0))./(fm^2+(f-f0).^2))./pi;  %--����������ɢ��
G_n=[G_f(nfft/2+1:end),G_f(1:nfft/2)];     
% G_n1=G_n(1:5);
rx=ifft(G_n);                          %�Թ����׽��и���Ҷ���任�����������rx(k)
rxk=rx(1:5);
[AR,E]=levinson(rxk,4);                                                      %---�õ�AR�˲�����ϵ��
  fc=37500+rand(1)*100;                                                  %---�źŵ�Ƶ��               
  %---����Ƶ��
   SNR=-40:0.1:0;
   length=size(SNR,2);
  alpha=0.01;
%-----�ڹ۲�ʱ����ͬ����£��Ƚϲ�ͬ�龯���ʺ�����ȶԼ���Ӱ��-----

 Time=0.0001;
  for i=1:length
        [ pd1(:,i), pa1(:,i)]=test_signal_2(fc,fs,SNR(i),Time,AR,alpha);              %---���ù��ܺ�����ʵ��
  end
%  
%   fs=20000000;
%   for i=1:length
%   [ pd2(:,i), pa2(:,i)]=test_signal_2(fc,fs,SNR(i),Time,AR,alpha);              %---���ù��ܺ�����ʵ��
%   end
%   fs=50000000;
%    for i=1:length
%   [ pd3(:,i), pa3(:,i)]=test_signal_2(fc,fs,SNR(i),Time,AR,alpha);              %---���ù��ܺ�����ʵ��
%    end
 
      figure(1)
    plot(SNR,pd1,'r')
    hold on
    grid on
 %     plot(SNR,pd2,'g')
 %     hold on
 %     plot(SNR,pd3,'b')
  title('��ɫ���������²�ͬ����Ƶ���������λ��Ƶ���ź�ʵ���龯����');
  xlabel('SNR(db)');
  ylabel('ϵͳ��ȷ����Pd');
%    legend('����Ƶ��Ϊ1MHZ','����Ƶ��20MHZ','����Ƶ��50MHZ');
  
     figure(2)
    plot(SNR,pa1,'r')
 %     hold on
 %     grid on
%      plot(SNR,pa2,'g')
%      hold on
 %     plot(SNR,pa3,'b')
  title('��ɫ���������²�ͬ����Ƶ���������λ��Ƶ���ź��龯������');
  xlabel('SNR(db)');
  ylabel('ϵͳ�龯����Pf');
 %   legend('����Ƶ��Ϊ1MHZ','����Ƶ��20MHZ','����Ƶ��50MHZ');
  
  
