%---�ǰ����������Ƶ�ʣ���λ�źŵļ��(�Ӻ���)
%---��ɫ�����µļ��

function [Pd,Pa]=test_signal_2(fc,fs,SNR,Time,AR,alpha)
   N=floor(fs*Time);                                %---�۲�ʱ����Ĳ�������
   f=fc/fs;                                         %---��һ���źŵ�Ƶ��
   A=1;                                             %---�̶��źŵķ�ֵ������ȵĸı��ͨ���ı������Ĺ�����ʵ��
   ph=2*pi*rand;
   S=A.*sin(2*pi*f*[1:N]+ph);
   Es = S*S';                                       %---�ź�����
   Ps=Es/N;
   Dn = Ps/power(10,SNR/10);                        %---��������
   Z0=1.6*sqrt(-Dn*Es*log(alpha));                   %---�о����ޣ����龯�����������������Z0ǰ��һϵ�����ı�Z0
   
%---------����������������5000��ʵ�������������---------��������������������
   Pd=0;                                            %---Pd���ڴ����ȷ������
   Pa=0;                                            %---Pa���ڴ���龯����
   jj=0;
 for k=1:500
   G_W_noise=wgn(1,N,1,'linear');                 %---������λ��˹������          
   noise=filter(1,AR,G_W_noise);                  %---������Ҫ�������״�ķǰ�����
  
   sum1=0;                                          %---��С�γ����ǰ��޸������Ĺ��ʣ�ʹ֮����SNR
   for i=1:N
       sum1=sum1+noise(i).^2;
   end
   sum1=sum1/N;
   para=Dn./sum1;
   noise=noise*sqrt(para);
   ph=2*pi*rand;
    if (rand>0.5)
             x=A.*sin(2*pi*f*[1:N]+ph)+noise;  
             ii=1;
             jj=jj+1;
      else
                 
             x=noise;                                        %���յ���ֻ������
            ii=0;
      end
 
   Ss=A.*sin(2*pi*f*[1:N]);                  
   Sc=A.*cos(2*pi*f*[1:N]);
   Gss(k)=x*Ss'; 
   Gcs(k)=x*Sc';
   Zs(k)=sqrt(Gss(k).^2+Gcs(k).^2);                
   if ((Zs(k)>Z0)&&(ii==1))
       Pd=Pd+1;
   end
   if ((Zs(k)>Z0)&&(ii==0))
  
       Pa=Pa+1;
   end
 end
  Pd=Pd/jj;
  Pa=Pa/(500-jj);

   


  
   


