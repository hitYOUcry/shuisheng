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
title('单通道和多通道下的虚警概率');
ylabel('虚警概率pf');
xlabel('信噪比snr');
legend('单通道','多通道');

figure(2)
plot(snr,pd1_est,'r');
hold on
grid on
plot(snr,pdarray_est);
title('单通道和多通道下的真确检测概率');
ylabel('检测概率pd');
xlabel('信噪比snr');
legend('单通道','多通道');


figure(3)
plot(snr,distance_est,'r');
hold on
grid on
axis([-21  0  1199  1201])
title('目标距离检测');
ylabel('目标距离');
xlabel('信噪比snr');

figure(4)
plot(snr,doa_est,'r');
hold on
grid on
title('目标方向检测');
ylabel('目标与阵列的角度');
xlabel('信噪比snr');

figure(5)
subplot(1,2,1);
plot([-90:90],Pbeam(6,:),'r');
hold on
stem(50,1,'b');
grid on
title('信噪比在-15db时，波束图');
ylabel('波束Pbeam');
xlabel('角度');

subplot(1,2,2);
plot([-90:90],Pbeam(16,:),'r');
hold on
grid on
stem(50,1,'b');
title('信噪比在-5db时，波束图');
ylabel('波束Pbeam');
xlabel('角度');









