%---�ǰ����������Ƶ�ʣ������λ�źŵļ��(�Ӻ���)
%---��ɫ�����µļ��

function  [Pd,Pa]=yubai(fs,SNR,Time,AR,alpha)
    N=floor(fs*Time);                                %---�۲�ʱ����Ĳ�������
    A=1;                                             %---�̶��źŵķ�ֵ������ȵĸı��ͨ���ı������Ĺ�����ʵ��
  
   
%---------����������������5000��ʵ�������������---------��������������������
   Pd=0;                                            %---Pd���ڴ����ȷ������
   Pa=0;                                            %---Pa���ڴ���龯����
    w=1:41;
    ph=2*pi*rand;
    S=A.*sin((2*pi*(220+w-1)/fs)'*[1:N]+ph); %ÿ�ж���һ��ʱ�����У���21��
   tem=S.*S; 
   for w=1:41
         Es(w,1)=sum(tem(w,:));
    end                                 %---�ź�����
    Ps=Es/N;
    Dn= Ps/power(10,SNR/10);                        %---�������ʣ���һ��������
    Z0=1.1*sqrt(-Dn.*Es*log(alpha));               %---�о����ޣ����龯�����������������Z0ǰ��һϵ�����ı�Z0����һ��������
   jj=0;
   for k=1:500
      f=220+floor(unifrnd(0,40));
       ph=2*pi*rand;
      s1=A.*sin(2*pi*f/fs*[1:N]+ph);
      Es1=s1*s1';
      Ps1=Es1/N;
      Dn1= Ps1/power(10,SNR/10);  
      G_W_noise=wgn(1,N,1,'linear');                 %---������λ��˹������    
    
       noise=filter(1,AR,G_W_noise);                  %---������Ҫ�������״�ķǰ�����
                                              %---��С�γ����ǰ��޸������Ĺ��ʣ�ʹ֮����SNR
        sum1=noise*noise';
       sum1=sum1/N;
       para=Dn1./sum1;
       noise=noise.*sqrt(para);
       
      
      if (rand>0.5)
             x=s1+noise;  
             ii=1;
             jj=jj+1;
      else
                 
             x=noise;                                        %���յ���ֻ������
            ii=0;
      end
      x=filter(AR,1,x) ;                               %---Ԥ�׻�
      
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
            Z_s(w)= Z_s(w)/Z0(w);               %�о�������Ӧ������ֵ��һ��   
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
 
  
     

