
%---��˹�������µļ��(�Ӻ���)

function Pd=test_signal_1(fc,ph,fs,SNR,Time,a)
   N=floor(fs*Time);                                 %---�۲�ʱ����Ĳ�������
   f=fc/fs;                                          %---��һ���źŵ�Ƶ��
   A=1;                                           
   S=A.*sin(2*pi*f*[1:N]+ph);
   Es = S*S';                                        %---�ź�����
   Ps=Es/N;
   Dn = Ps/power(10,SNR/10);                         %---��������
   Z0=1.1*sqrt(-2*Dn*Es*log(a));                     %---�о�����
 

%---------��5000��ʵ�������������---------
   Pd=0;                                             %---Pd���ڴ����ȷ������
   Pa=0;                                             %---Pa���ڴ���龯����
 for k=1:5000
   noise=wgn(1,N,Dn,'linear');                       %---��������ΪDn(watt)�ĸ�˹������          
   x=A.*sin(2*pi*f*[1:N]+ph)+noise;                  %---���յ������źŵ��������
   Ss=A.*sin(2*pi*f*[1:N]);                          %---ȷ���źŵ�Ƶ��   
   Sc=A.*cos(2*pi*f*[1:N]);
   Gss(k)=x*Ss'; 
   Gcs(k)=x*Sc';
   Zs(k)=sqrt(Gss(k).^2+Gcs(k).^2);                  %---�źŵ�Z�任
   if Zs(k)>Z0
       Pd=Pd+1;
   end

   Gsn(k)=noise*Ss';                                %---���յ����źŵ��������
   Gsc(k)=noise*Sc';
   Zn(k)=sqrt(Gsn(k).^2+Gsc(k).^2);
   if Zn(k)>Z0
       Pa=Pa+1;
   end
 end
  Pd=Pd/5000;
  Pa=Pa/5000;
  
  
  
   


