%����Ϊ��˹���������ź�Ƶ��ȷ��
function Pd=WGN_phasetest(fs,fc,a,SNR,Tlen,N_Exp)

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
