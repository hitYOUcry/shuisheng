%背景为高斯白噪声且信号频率确定
function Pd=WGN_phasetest(fs,fc,a,SNR,Tlen,N_Exp)

N = fs*Tlen;                        %采样点数
fc_norm = fc/fs;                    %归一化的目标信号频率
A = 1;                              %信号幅度
Ps = A^2/2;                         %信号功率
Es = Ps*N;                          %信号能量  P55 4－33
Dn = Ps/power(10,SNR/10);           %噪声功率，即其方差

gs = A*sin(2*pi*fc_norm*[0:N-1])';  %相关器参考信号
gc = A*cos(2*pi*fc_norm*[0:N-1])';  %相关器参考信号

Zt = sqrt(-2*Dn*Es*log(a));         %判决门限   P58 4－46
Pd1 = marcumq(sqrt(Es/Dn),sqrt(-2*log(a)));       %检测概率理论值   
%Q_m(a,b) = 1/a^(m-1) * integral from b Tlen inf of[x^m * exp(-(x^2+a^2)/2) * I_(m-1)(ax)] dx.
%a=d m=1 b=z0 z0见P58 4－46


Pd = 0;
Pa = 0;
for(k=1:N_Exp)                                      %5000次实验
  s = A*sin(2*pi*fc_norm*[0:N-1]+2*pi*rand)';       %目标信号  rand产生一个随机数 代表随机相位
  n = wgn(N,1,Dn,'linear');         %WGN噪声
  Gs = gs'*(s+n);                   %H1假设时的相关器输出
  Gc = gc'*(s+n);                   %H1假设时的相关器输出
  Z1(k) = sqrt(Gs^2+Gc^2);          %H1假设时的Z计算器输出
  if(Z1(k) > Zt)                    %判决
    Pd = Pd+1;
  end
  Gs = gs'*n;                       %H0假设时的相关器输出
  Gc = gc'*n;                       %H0假设时的相关器输出
  Z0(k) = sqrt(Gs^2+Gc^2);          %H0假设时的Z计算器输出
  if(Z0(k) > Zt)                    %判决
    Pa = Pa+1;
  end
end
Pd = Pd/N_Exp;                       %蒙特卡罗实验得出的实际检测概率
Pa = Pa/N_Exp;

%-----------------------------------------------------------------------%
