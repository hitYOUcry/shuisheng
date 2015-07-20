%����Ϊ��˹���������ź�Ƶ��ȷ��

%-----------------------------------------------------------------------%
%��ʼ��
fs = 4000;                          %����Ƶ��(Hz)
fc = 250;                           %�ź�Ƶ��(Hz)
a = 0.01;                           %�龯��������
SNR=5;                             %�����
Tlen=250/1000;                       %�۲�ʱ��T
N_Exp=5000;                         %ʵ�����
%-----------------------------------------------------------------------%


N = fs*Tlen;                        %��������
fc_norm = fc/fs;                    %��һ����Ŀ���ź�Ƶ��
A = 1;                              %�źŷ���
Ps = A^2/2;                         %�źŹ���
Es = Ps*N;                          %�ź�����  P55 4��33
Dn = Ps/power(10,SNR/10);           %�������ʣ����䷽��

gs = A*sin(2*pi*fc_norm*[0:N-1])';  %������ο��ź�
gc = A*cos(2*pi*fc_norm*[0:N-1])';  %������ο��ź�

Zt = sqrt(-2*Dn*Es*log(a));         %�о�����   P58 4��46
Pd1 = marcumq(sqrt(Es/Dn),sqrt(-2*log(a)));       %����������ֵ   
%Q_m(a,b) = 1/a^(m-1) * integral from b Tlen inf of[x^m * exp(-(x^2+a^2)/2) * I_(m-1)(ax)] dx.
%a=d m=1 b=z0 z0��P58 4��46


Pd = 0;
Pa = 0;
for(k=1:N_Exp)                                      %5000��ʵ��
  s = A*sin(2*pi*fc_norm*[0:N-1]+2*pi*rand)';       %Ŀ���ź�  rand����һ������� ���������λ
  n = wgn(N,1,Dn,'linear');         %WGN����
  Gs = gs'*(s+n);                   %H1����ʱ����������
  Gc = gc'*(s+n);                   %H1����ʱ����������
  Z1(k) = sqrt(Gs^2+Gc^2);          %H1����ʱ��Z���������
  if(Z1(k) > Zt)                    %�о�
    Pd = Pd+1;
  end
  Gs = gs'*n;                       %H0����ʱ����������
  Gc = gc'*n;                       %H0����ʱ����������
  Z0(k) = sqrt(Gs^2+Gc^2);          %H0����ʱ��Z���������
  if(Z0(k) > Zt)                    %�о�
    Pa = Pa+1;
  end
end
Pd = Pd/N_Exp;                       %���ؿ���ʵ��ó���ʵ�ʼ�����
Pa = Pa/N_Exp;

%-----------------------------------------------------------------------%
%���ʵ����
figure;                             
subplot(2,1,1);
plot(Z1);
hold on;
plot(Zt*ones(1,N_Exp),'r');
title({'��˹�������������λ�źŵļ��(5000��ʵ��)',...
      ['����Ƶ��fs(Hz) = ',num2str(fs),', �ź�Ƶ��fc(Hz) = ',num2str(fc),...
       ', �۲�ʱ��Tlen(ms) = ',num2str(Tlen*1000), ', SNR(dB) = ',num2str(SNR)]});
text(2500,max(Z1)-200,{['���ۼ�����Pd1 = ',num2str(Pd1)],['ʵ�������Pd = ',num2str(Pd)],...
    ['�о�����Zt = ',num2str(Zt)]},'FontSize',8);
legend('�����Z1','�о�����Zt');
xlabel('ʵ�����k','FontSize',10);
ylabel('H1����ʱ�ļ����Z1','FontSize',10);
                           
subplot(2,1,2);
plot(Z0);
hold on;
plot(Zt*ones(1,5000),'r');
text(2500,max(Z0)-2,{['�龯��������a = ',num2str(a)],['ʵ���龯����Pa = ',num2str(Pa)],...
     ['�о�����Zt = ',num2str(Zt)]},'FontSize',8);
legend('�����Z0','�о�����Zt');
xlabel('ʵ�����k','FontSize',10);
ylabel('H0����ʱ�ļ����Z0','FontSize',10);