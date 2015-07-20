%背景为高斯白噪声且信号频率确定

%-----------------------------------------------------------------------%
%初始化
fs = 4000;                          %采样频率(Hz)
fc = 250;                           %信号频率(Hz)
a = 0.01;                           %虚警概率上限
SNR=5;                             %信噪比
Tlen=250/1000;                       %观测时间T
N_Exp=5000;                         %实验次数
%-----------------------------------------------------------------------%


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
%输出实验结果
figure;                             
subplot(2,1,1);
plot(Z1);
hold on;
plot(Zt*ones(1,N_Exp),'r');
title({'高斯白噪声中随机相位信号的检测(5000次实验)',...
      ['采样频率fs(Hz) = ',num2str(fs),', 信号频率fc(Hz) = ',num2str(fc),...
       ', 观测时间Tlen(ms) = ',num2str(Tlen*1000), ', SNR(dB) = ',num2str(SNR)]});
text(2500,max(Z1)-200,{['理论检测概率Pd1 = ',num2str(Pd1)],['实验检测概率Pd = ',num2str(Pd)],...
    ['判决门限Zt = ',num2str(Zt)]},'FontSize',8);
legend('检测量Z1','判决门限Zt');
xlabel('实验次数k','FontSize',10);
ylabel('H1假设时的检测量Z1','FontSize',10);
                           
subplot(2,1,2);
plot(Z0);
hold on;
plot(Zt*ones(1,5000),'r');
text(2500,max(Z0)-2,{['虚警概率上限a = ',num2str(a)],['实验虚警概率Pa = ',num2str(Pa)],...
     ['判决门限Zt = ',num2str(Zt)]},'FontSize',8);
legend('检测量Z0','判决门限Zt');
xlabel('实验次数k','FontSize',10);
ylabel('H0假设时的检测量Z0','FontSize',10);