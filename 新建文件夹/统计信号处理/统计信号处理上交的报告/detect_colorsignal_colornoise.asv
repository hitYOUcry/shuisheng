function [pf_test,pd_test]=detect_colorsignal_colornoise(T,snr,pf,G0)
%%%detect_colorsignal_colornoise: 高斯非白噪声背景下检测高斯信号---Monte Carlo 仿真
%%%input                    T-----------观测时间
%%%                         snr---------信噪比(信号平均功率与噪声平均功率之比)，用数值表示10^(snr/10)
%%%                         pf----------限定的虚警概率  pf<=0.01
%%%                         G0----------仿真门限
%%%output                  
%%%                         pf_test-----仿真实验得到的虚警概率
%%%                         pd_test-----仿真实验得到的正确检测概率
%%%by tanjunhong,on 7.2,2011
%%%modified on 


% T=0.01,snr=0,pf=0.01    %test this function

%%%----信号参数设置----------------------------------------
fs=4000;                  %对接收数据的采样频率
ts=1/fs;                  %采样时间间隔
N=T/ts;                   %观测数据个数(时域、频域个数一样)
t=(1:N)'*ts;              %采样时刻的列矢量
df=1/T;                   %采样频域间隔
%%%----色噪声、信号功率谱参数设置---------------------------------
k_n=0.1;fm_n=200;f0_n=260; %色噪声功率谱G(f)参数选择
k_s=1;fm_s=400;f0_s=360;   %色信号功率谱G(f)参数选择
snr=10^(snr/10);           %信噪比

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%----由非白噪声、非白高斯信号的功率谱确定厄卡特滤波器-----------------------------
f=-fs/2:df:fs/2-df;        %共N个频率采样点
f=f';                      %列矢量

PNW0=1/pi.*((fm_n+k_n*(f+f0_n))./(fm_n^2+(f+f0_n).^2)+(fm_n-k_n*(f-f0_n))./(fm_n^2+(f-f0_n).^2))*fs;   % 给定的功率谱G(f)：[-pi,pi]
PNW=[PNW0(N/2+1:end);PNW0(1:N/2)];                                                   %[0,2*pi]
PSW0=1/pi.*((fm_s+k_s*(f+f0_s))./(fm_s^2+(f+f0_s).^2)+(fm_s-k_s*(f-f0_s))./(fm_s^2+(f-f0_s).^2))*fs;   % 给定的功率谱G(f)：[-pi,pi]
PSW=[PSW0(N/2+1:end);PSW0(1:N/2)];                                               %[0,2*pi]

HW0=1./PNW-1./(PNW+PSW);
HW=sqrt(HW0);                      %厄卡特滤波器

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%  确定仿真门限     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G0=13;                             %观测时间0.01
% % % G0=0.9e4;                          %观测时间0.1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   Monte Carlo  仿真   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pf_test=0;
pd_test=0;
N_exp=5000;                          %仿真次数
for ii=1:N_exp
    %%%%%%       仿真非白噪声信号     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wn=color_produce(k_n,fm_n,f0_n,N);%列矢量
    Pn=wn'*wn/N;                      %噪声功率
    %%%%%%  仿真非白高斯信号  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s0=color_produce(k_s,fm_s,f0_s,N);%列矢量
    Ps=s0'*s0/N;
    Ds=snr*Pn;                        %需要的信号功率
    s=s0*sqrt(Ds/Ps);                 %符合信噪比的信号
    
    
    %%% H1假设下
    x1=s+wn;
    XW1=fft(x1);
    YW1=XW1.*HW;
    G_1=YW1'*YW1*T;                     %检测统计量
% %     y1=ifft(YW1);
% %     G_1=y1'*y1*T;
    if G_1>G0
        pd_test=pd_test+1;
    end
    
    %%%% H0假设下
    x0=wn;
    XW0=fft(x0);
    YW0=XW0.*HW;
    G_0=YW0'*YW0*T;                     %检测统计量
% %     y1=ifft(YW1);
% %     G_1=y1'*y1*T;
    if G_0>G0
        pf_test=pf_test+1;
    end
    
end
pf_test=pf_test/N_exp;
pd_test=pd_test/N_exp;

% pd_test,pf_test          % to test this function

% %%%%%%%%%%%try to find the threhold G0:13 （T=0.01）, 0.9e4 (T=0.1)
% figure
% plot(1:5000,G_1,'r*');hold on;grid on
% plot(1:5000,G_0,'b*');









