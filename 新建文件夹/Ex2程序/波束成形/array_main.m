clc
clear all;
snr=-20:0;
T=0.01;
pf=0.01;
for jj=1:21
    [pf1_est(jj),pfarray_est(jj),pd1_est(jj),pdarray_est(jj),distance_est(jj),doa_est(jj),Pbeam(jj,:)]=detect_array_know_WGN(T,snr(jj),pf);
end
figure(1)
plot(snr,pf1_est,'r');
hold on
grid on
axis([-21 0  0  0.1])
plot(snr,pfarray_est);
title('��ͨ���Ͷ�ͨ���µ��龯����');
ylabel('�龯����pf');
xlabel('�����snr');
legend('��ͨ��','��ͨ��');

figure(2)
plot(snr,pd1_est,'r');
hold on
grid on
plot(snr,pdarray_est);
title('��ͨ���Ͷ�ͨ���µ���ȷ������');
ylabel('������pd');
xlabel('�����snr');
legend('��ͨ��','��ͨ��');


figure(3)
plot(snr,distance_est,'r');
hold on
grid on
axis([-21  0  1199  1201])
title('Ŀ�������');
ylabel('Ŀ�����');
xlabel('�����snr');

figure(4)
plot(snr,doa_est,'r');
hold on
grid on
title('Ŀ�귽����');
ylabel('Ŀ�������еĽǶ�');
xlabel('�����snr');

figure(5)
subplot(1,2,1);
plot([-90:90],Pbeam(6,:),'r');
hold on
stem(50,1,'b');
grid on
title('�������-15dbʱ������ͼ');
ylabel('����Pbeam');
xlabel('�Ƕ�');

subplot(1,2,2);
plot([-90:90],Pbeam(16,:),'r');
hold on
grid on
stem(50,1,'b');
title('�������-5dbʱ������ͼ');
ylabel('����Pbeam');
xlabel('�Ƕ�');









