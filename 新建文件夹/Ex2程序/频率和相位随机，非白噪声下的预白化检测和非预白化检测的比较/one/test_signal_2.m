%---�ǰ����������Ƶ�ʣ���λ�źŵļ��(�Ӻ���)
%---��ɫ�����µļ��

function  [Pd,Pa]=test_signal_2(fc,fs,SNR,Time,AR,alpha)
   N=floor(fs*Time);                                %---�۲�ʱ����Ĳ�������
   A=1;                                             %---�̶��źŵķ�ֵ������ȵĸı��ͨ���ı������Ĺ�����ʵ��
   f=fc/fs;
   ph=2*pi*rand;
   S=A.*sin(2*pi*(240/fs)*[1:N]+ph);
   Es = S*S';                                       %---�ź�����
   Ps=Es/N;
   Dn = Ps/power(10,SNR/10);                        %---��������
   Z0=1.85* sqrt(-Dn*Es*log(alpha));                   %---�о����ޣ����龯�����������������Z0ǰ��һϵ�����ı�Z0
   
%---------����������������5000��ʵ�������������---------��������������������
   Pd=0;                                            %---Pd���ڴ����ȷ������
   Pa=0;                                            %---Pa���ڴ���龯����
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
   x=A.*sin(2*pi*f*[1:N]+ph)+noise;                 %---���յ������źŵ��������
   for w=0:10
    
        Ss=A*sin(2*pi*((220+w*4)/fs)*[1:N]);                  
        Sc=A*cos(2*pi*((220+w*4)/fs)*[1:N]);
        Gss(k,w+1)=x*Ss'; 
        Gcs(k,w+1)=x*Sc';
        Z_s(k,w+1)=sqrt(Gss(k,w+1).^2+Gcs(k,w+1).^2);                  %---�źŵ�Z�任
        
       Gsn(k,w+1)=noise*Ss';                                 %---���յ������źŵ��������
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
 
  
   


