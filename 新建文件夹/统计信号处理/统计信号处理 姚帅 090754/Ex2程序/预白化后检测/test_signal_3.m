%---�ǰ����������Ƶ�ʣ���λ�źŵļ��(�Ӻ���)

function  Pd=test_signal_3(fc,ph,fs,SNR,Time,AR,a)
   N=floor(fs*Time);                                %---�۲�ʱ����Ĳ�������
   f=fc/fs;                                         %---��һ���źŵ�Ƶ��
   A=1;                                             %---�̶��źŵķ�ֵ������ȵĸı��ͨ���ı������Ĺ�����ʵ��
   S=A.*sin(2*pi*f*[1:N]+ph);
   Es = S*S';                                       %---�ź�����
   Ps=Es/N;
   Dn = Ps/power(10,SNR/10);                        %---��������
   Z0= 0.5*sqrt(-2*Dn*Es*log(a));                   %---�о����ޣ����龯�����������������Z0ǰ��һϵ�����ı�Z0
   
%---------����������������5000��ʵ�������������---------��������������������
   Pd=0;                                            %---Pd���ڴ����ȷ������
   Pa=0;                                            %---Pa���ڴ���龯����
 for k=1:5000
   G_W_noise=wgn(1,N,1,'linear');                  %---������λ��˹������          
   noise=filter(1,AR,G_W_noise);                   %---������Ҫ�������״�ķǰ�����
  
   sum1=0;                                          %---��С�γ����ǰ��޸������Ĺ��ʣ�ʹ֮����SNR
   for i=1:N
       sum1=sum1+noise(i).^2;
   end
   sum1=sum1/N;
   para=Dn./sum1;
   noise=noise*sqrt(para);
  
   x=A.*sin(2*pi*f*[1:N]+ph)+noise;                 %---���յ������źŵ��������
   x=filter(AR,1,x);                                %---Ԥ�׻�
   Ss=A.*sin(2*pi*f*[1:N]);                  
   Sc=A.*cos(2*pi*f*[1:N]);
   Gss(k)=x*Ss'; 
   Gcs(k)=x*Sc';
   Zs(k)=sqrt(Gss(k).^2+Gcs(k).^2);               
   if Zs(k)>Z0
       Pd=Pd+1;
   end
   
   noise=filter(AR,1,noise);
   Gsn(k)=noise*Ss';                                 %---���յ������źŵ��������
   Gsc(k)=noise*Sc';
   Zn(k)=sqrt(Gsn(k).^2+Gsc(k).^2);
   if Zn(k)>Z0
       Pa=Pa+1;
   end
 end
  Pd=Pd/5000;
  Pa=Pa/5000;
