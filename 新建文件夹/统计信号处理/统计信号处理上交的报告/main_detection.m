function main_detection()
%%%Main_detection:   plot the ROC curve
%%%by tanjunhong,on 6.30,2011
%%%modified on 7.1,2011
%%%modified on

clc
clear all
close all


snr=-15:1:5;           %�����ʸ��--ROC���ߺ�����
T=[0.01,0.1];          %�۲�ʱ��10ms,100ms
pf=0.01;               %�޶��龯����

%%    ȷ֪�źŵļ�⣨WGN��

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%    ȷ֪�źŵļ�⣨WGN��������ROC curve   %%%%%%%%%%%%%%%%%%%%%%
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
% legend('�۲�ʱ��Ϊ10ms','�۲�ʱ��Ϊ100ms');
% title('ȷ֪�źŵļ��     ���ջ��Ĺ�����������');

%%    �����λ�źŵļ�⣨WGN��

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%    �����λ�źŵļ�⣨WGN��������ROC curve   %%%%%%%%%%%%%%%%%%
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
% legend('�۲�ʱ��Ϊ10ms','�۲�ʱ��Ϊ100ms');
% title('�����λ�źŵļ��     ���ջ��Ĺ�����������');

%%    �����λ��Ƶ���źŵļ�⣨WGN��

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%    �����λ��Ƶ���źŵļ�⣨WGN��������ROC curve   %%%%%%%%%%%%
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
% legend('�۲�ʱ��Ϊ10ms','�۲�ʱ��Ϊ100ms');
% title('�����λ��Ƶ���źŵļ��     ���ջ��Ĺ�����������');
% 
% figure
% plot(snr,pf_test1,'b-*');hold on
% plot(snr,pf_test2,'r-*');
% plot(snr,pf,'k-*');grid on
% xlabel('SNR(dB)');ylabel('�龯����  Pf');
% legend('�۲�ʱ��Ϊ10ms','�۲�ʱ��Ϊ100ms','�㶨���龯����');
% title('�����λ��Ƶ���źŵļ��     �龯����');

%%    �����λ��Ƶ���źŵļ�⣨��ɫ������

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% �����λ��Ƶ���źŵļ�⣨��ɫ������������ROC curve���Ƚ�Ԥ�״������Ԥ�״���   %%%%%%%%
% %%%%%%%%%%%%      �۲�ʱ��10ms���龯����Pf=0.01        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %��Ԥ�״���
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for ii=1:length(snr)
%     [pf_test1(ii),pd_test1(ii),pd1(ii)]=detect_pf_unknow_colornoise(T(1),snr(ii),pf);
% end
% %Ԥ�״���
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
% legend('��Ԥ�״���','Ԥ�״���1','Ԥ�״���2');
% title('�����λ��Ƶ���źŵļ�⣨ɫ������   ���ջ��Ĺ�����������');
% figure
% plot(snr,pf_test1,'b-*');hold on;grid on
% plot(snr,pf_test2,'r-*');
% plot(snr,pf_test3,'g-*');
% pf_constant=pf*ones(1,length(snr));
% plot(snr,pf_constant,'k-*');grid on
% xlabel('SNR(dB)');ylabel('�龯����  Pf');
% legend('��Ԥ�״���','Ԥ�״���1','Ԥ�״���2','�㶨���龯����');
% title('�����λ��Ƶ���źŵļ�⣨ɫ������   �龯����');


%%    �ǰ������·ǰ׸�˹�źŵļ������ʱ��ó�20����

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%    �ǰ������·ǰ׸�˹�źŵļ�⡪����ROC curve      %%%%%%%%%%%%
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
% legend('�۲�ʱ��Ϊ10ms','�۲�ʱ��Ϊ100ms');
% title('��˹���������¸�˹�źŵļ��     ���ջ��Ĺ�����������');
% figure
% plot(snr,pf_test1,'b-*');hold on
% plot(snr,pf_test2,'r-*');
% pf_constant=pf*ones(1,length(snr));
% plot(snr,pf_constant,'k-*');grid on
% xlabel('SNR(dB)');ylabel('�龯����  Pf');
% legend('�۲�ʱ��Ϊ10ms','�۲�ʱ��Ϊ100ms','�㶨���龯����');
% title('��˹���������¸�˹�źŵļ��    �龯����');

%%    ����������̽��ϵͳ���棨ȷ֪�źŵļ�⣩Elapsed time is 4067.452734 seconds.(67����)


















