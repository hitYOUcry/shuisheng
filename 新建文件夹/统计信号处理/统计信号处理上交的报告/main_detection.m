function main_detection()
%%%Main_detection:   plot the ROC curve
%%%by tanjunhong,on 6.30,2011
%%%modified on 7.1,2011
%%%modified on

clc
clear all
close all


snr=-15:1:5;           %信噪比矢量--ROC曲线横坐标
T=[0.01,0.1];          %观测时间10ms,100ms
pf=0.01;               %限定虚警概率

%%    确知信号的检测（WGN）

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%    确知信号的检测（WGN）———ROC curve   %%%%%%%%%%%%%%%%%%%%%%
% for ii=1:length(snr)
%     [pf_test1(ii),pd_test1(ii),pd1(ii)]=detect_known_WGN(T(1),snr(ii),pf);
% end
% for ii=1:length(snr)
%     [pf_test2(ii),pd_test2(ii),pd2(ii)]=detect_known_WGN(T(2),snr(ii),pf);
% end
% 
% figure
% plot(snr,pd_test1,'b');hold on
% plot(snr,pd_test2,'r');grid on
% xlabel('SNR(dB)');ylabel('detection probability  Pd');
% legend('观测时间为10ms','观测时间为100ms');
% title('确知信号的检测     接收机的工作特性曲线');

%%    随机相位信号的检测（WGN）

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%    随机相位信号的检测（WGN）———ROC curve   %%%%%%%%%%%%%%%%%%
% for ii=1:length(snr)
%     [pf_test1(ii),pd_test1(ii),pd1(ii)]=detect_p_unknow_WGN(T(1),snr(ii),pf);
% end
% for ii=1:length(snr)
%     [pf_test2(ii),pd_test2(ii),pd2(ii)]=detect_p_unknow_WGN(T(2),snr(ii),pf);
% end
% 
% figure
% plot(snr,pd_test1,'b');hold on
% plot(snr,pd_test2,'r');grid on
% xlabel('SNR(dB)');ylabel('detection probability  Pd');
% legend('观测时间为10ms','观测时间为100ms');
% title('随机相位信号的检测     接收机的工作特性曲线');

%%    随机相位、频率信号的检测（WGN）

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%    随机相位、频率信号的检测（WGN）———ROC curve   %%%%%%%%%%%%
% for ii=1:length(snr)
%     [pf_test1(ii),pd_test1(ii),pd1(ii)]=detect_pf_unknow_WGN(T(1),snr(ii),pf);
% end
% for ii=1:length(snr)
%     [pf_test2(ii),pd_test2(ii),pd2(ii)]=detect_pf_unknow_WGN(T(2),snr(ii),pf);
% end
% 
% figure
% plot(snr,pd_test1,'b');hold on
% plot(snr,pd_test2,'r');grid on
% xlabel('SNR(dB)');ylabel('detection probability  Pd');
% legend('观测时间为10ms','观测时间为100ms');
% title('随机相位、频率信号的检测     接收机的工作特性曲线');
% 
% figure
% plot(snr,pf_test1,'b-*');hold on
% plot(snr,pf_test2,'r-*');
% plot(snr,pf,'k-*');grid on
% xlabel('SNR(dB)');ylabel('虚警概率  Pf');
% legend('观测时间为10ms','观测时间为100ms','恒定的虚警概率');
% title('随机相位、频率信号的检测     虚警概率');

%%    随机相位、频率信号的检测（有色噪声）

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 随机相位、频率信号的检测（有色噪声）———ROC curve（比较预白处理与非预白处理）   %%%%%%%%
% %%%%%%%%%%%%      观测时间10ms，虚警概率Pf=0.01        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %非预白处理
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for ii=1:length(snr)
%     [pf_test1(ii),pd_test1(ii),pd1(ii)]=detect_pf_unknow_colornoise(T(1),snr(ii),pf);
% end
% %预白处理
% %% method 1----invR=D'*D;
% for ii=1:length(snr)
%     [pf_test2(ii),pd_test2(ii),pd2(ii)]=detect_pf_unknow_colornoise_yubai(T(1),snr(ii),pf);
% end
% %% method 2---filter(AR,b,wn)  
% for ii=1:length(snr)
%     [pf_test3(ii),pd_test3(ii),pd3(ii)]=detect_pf_unknow_colornoise_yubai_filter(T(1),snr(ii),pf);
% end
% figure(1)
% plot(snr,pd_test1,'b');hold on
% plot(snr,pd_test2,'r');
% plot(snr,pd_test3,'g');grid on
% xlabel('SNR(dB)');ylabel('detection probability  Pd');
% legend('非预白处理','预白处理1','预白处理2');
% title('随机相位、频率信号的检测（色噪声）   接收机的工作特性曲线');
% figure
% plot(snr,pf_test1,'b-*');hold on;grid on
% plot(snr,pf_test2,'r-*');
% plot(snr,pf_test3,'g-*');
% pf_constant=pf*ones(1,length(snr));
% plot(snr,pf_constant,'k-*');grid on
% xlabel('SNR(dB)');ylabel('虚警概率  Pf');
% legend('非预白处理','预白处理1','预白处理2','恒定的虚警概率');
% title('随机相位、频率信号的检测（色噪声）   虚警概率');


%%    非白噪声下非白高斯信号的检测运行时间好长20分钟

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%    非白噪声下非白高斯信号的检测———ROC curve      %%%%%%%%%%%%
% for ii=1:length(snr)
%     [pf_test1(ii),pd_test1(ii)]=detect_colorsignal_colornoise(T(1),snr(ii),pf,14);
% end
% for ii=1:length(snr)
%     [pf_test2(ii),pd_test2(ii)]=detect_colorsignal_colornoise(T(2),snr(ii),pf,9.5e3);
% end
% 
% figure
% plot(snr,pd_test1,'b');hold on
% plot(snr,pd_test2,'r');grid on
% xlabel('SNR(dB)');ylabel('detection probability  Pd');
% legend('观测时间为10ms','观测时间为100ms');
% title('高斯噪声背景下高斯信号的检测     接收机的工作特性曲线');
% figure
% plot(snr,pf_test1,'b-*');hold on
% plot(snr,pf_test2,'r-*');
% pf_constant=pf*ones(1,length(snr));
% plot(snr,pf_constant,'k-*');grid on
% xlabel('SNR(dB)');ylabel('虚警概率  Pf');
% legend('观测时间为10ms','观测时间为100ms','恒定的虚警概率');
% title('高斯噪声背景下高斯信号的检测    虚警概率');

%%    线阵列主动探测系统仿真（确知信号的检测）Elapsed time is 4067.452734 seconds.(67分钟)


















