function [pf1_est,pfarray_est,pd1_est,pdarray_est,distance_est,doa_est,Pbeam]=detect_array_know_WGN(T,snr,pf)
%%%%detect_array_know_WGN:线阵列主动探测系统-多元阵确知信号的检测
%%%%input：               T------------观测时间
%%%%                      snr----------接收机端的信噪比（处理前）
%%%%                      pf-----------限定的虚警概率 
%%%%output:               pd1_est------单通道确知信号的正确检测概率
%%%%                      pdarray_est--多通道确知信号的正确检测概率
%%%%                      distance_est--估计的目标距离
%%%%                                    需调用函数matchedfilter计算第一个基元接收信号相对发射信号的延迟时间从而计算出目标距离
%%%%                      doa_est-------估计的目标方位
%%%                       Pbeam---------空间谱（波束图）
%%%%by tanjunhong,7.5,2011（the whole day, from moring to 20:30,the run this function for one more hour）
%%%%modified on 

% T=0.01;snr=-15;pf=0.01;  % test this function   M=10,50

%%%%%%%%%%%%%%%%参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%信号参数
fs=4000;                  %对接收数据的采样频率
ts=1/fs;                  %采样时间间隔
N=T/ts;                   %观测数据个数
t=(1:N)*ts;               %采样时刻的行矢量
c=1500;                   %声速
A=1;                      %信号幅度
fc=250;                   %信号频率
lamda=c/fc;               %信号波长
phi=pi/4;                 %信号相位
s=A*sin(2*pi*fc*t+phi);   %行信号矢量------仿真确知信号

%%噪声参数
snr=10^(snr/10);          %信噪比
Ps=A^2/2;                 %信号的平均功率
Pn=Ps/snr;                %噪声的平均功率
sigma=sqrt(Pn);           %噪声的标准差（均值为0）

%%目标参数
distance=1200;            %目标与第一个基元距离
doa=50;                   %目标方位（与阵列法向的夹角）-90:90度
delay0=2*distance/c;      %第一个基元经过延迟tau0接收到主动声呐发射信号
delay0_N=round(delay0/ts);
%%阵列参数
M=40;                     %基元数  M=10
d=lamda/2;                %相邻基元间隔
tau0=(0:M-1)'*d*sin(doa*pi/180)/c;  %各个基元的延迟时间（相对参考基元1）
atheta0=exp(-1i*2*fc*tau0);         %导向矢量
theta=-90:90;             %阵列扫描角度范围
theta_N=length(theta);    %扫描角度数目


%%%-理论的正确检测概率------------------------------------------------------
d1=snr*N;                 %单通道
darray=snr*N*M;           %多通道
u0=qfuncinv(pf);          %由虚警概率pf反求u0  
pd1=qfunc(u0-sqrt(d1));        %单通道：理论的正确检测概率值
pdarray=qfunc(u0-sqrt(darray));%多通道：理论的正确检测概率


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  确定仿真门限     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G01=u0*sqrt(d1);           %单通道
G0array=u0*sqrt(darray);   %多通道


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   Monte Carlo  仿真   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pbeam_temp=0;
pd1_est=0;
pdarray_est=0;
pf1_est=0;
pfarray_est=0;
doa_est=0;
distance_est=0;
N_exp=50;                %仿真的实验次数

for ii=1:N_exp
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % %%%%%%%%%%%H1假设下（有信号）%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r=atheta0*s+sigma*randn(M,N);
    %%%% 单通道检测
%     G_1_1=r(1,:)*s'/Pn;    
    G_1_1=(s+sigma*randn(1,N))*s'/Pn;                        %单通道检测统计量
    if G_1_1>G01                               %比较判决
        pd1_est=pd1_est+1;
    end
    %%%%阵列信号检测,在有信号时，进行目标距离估计distance_est和目标方位估计doa_est
    for k=1:theta_N
        tau=(0:M-1)'*d*sin(theta(k)*pi/180)/c; %各个基元的延迟时间（相对参考基元1）
        atheta=exp(-1i*2*fc*tau);              %导向矢量M*1
        beam_out(k,:)=(atheta')*r;         %波束形成，对各个基元的接收信号进行延迟相加
        Pbeam(k)=beam_out(k,:)*beam_out(k,:)'; %空间谱，波束图 
    end
    Pbeam_temp=Pbeam+Pbeam_temp;
    [Pmax,theta_index]=max(Pbeam);
    doa_est=doa_est+theta(theta_index);        %doa估计
    
    G_array_1=beam_out(theta_index,:)*s'/Pn;   %阵列检测统计量
    if G_array_1>G0array                       %比较判决
        pdarray_est=pdarray_est+1;              
    end
    
                                               %目标距离估计
    r1=[zeros(1,delay0_N),r(1,:)];
    smf1=matchedfilter(r1,s,fs);
    [smfmax,delay_index]=max(smf1);
    distance_est=distance_est+delay_index*ts*c/2;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
% % %%%%%%%%%%%H0假设下（无信号）%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r=sigma*randn(M,N);
    %%%% 单通道检测
    G_1_0=r(1,:)*s'/Pn;                        %单通道检测统计量
    if G_1_0>G01                               %比较判决
        pf1_est=pf1_est+1;
    end
    %%%%阵列信号检测,在有信号时，进行目标距离估计distance_est和目标方位估计doa_est
    for k=1:theta_N
        tau=(0:M-1)'*d*sin(theta(k)*pi/180)/c;  %各个基元的延迟时间（相对参考基元1）
        atheta=exp(-1i*2*fc*tau);               %导向矢量M*1
        beam_out0(k,:)=(atheta')*r;              %波束形成，对各个基元的接收信号进行延迟相加
        Pbeam0(k)=beam_out0(k,:)*beam_out0(k,:)';  %空间谱，波束图 
    end
    Pbeam_temp0=Pbeam0+Pbeam_temp;
    [Pmax0,theta_index0]=max(Pbeam0);
    G_array_0=beam_out0(theta_index0,:)*s'/Pn;   %阵列检测统计量
    if G_array_0>G0array                       %比较判决
        pfarray_est=pfarray_est+1;              
    end 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
    
end
Pbeam=Pbeam_temp/N_exp;
Pbeam=Pbeam/max(Pbeam);
pd1_est=pd1_est/N_exp;
pdarray_est=pdarray_est/N_exp;
pf1_est=pf1_est/N_exp;
pfarray_est=pfarray_est/N_exp;
doa_est=doa_est/N_exp;
distance_est=distance_est/N_exp;



% pd1,pd1_est,pdarray,pdarray_est,doa_est,distance_est,pf1_est,pfarray_est



   


