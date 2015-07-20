function [pf1_est,pfarray_est,pd1_est,pdarray_est,distance_est,doa_est,Pbeam]=detect_array_know_WGN(T,snr,pf)
%%%%detect_array_know_WGN:����������̽��ϵͳ-��Ԫ��ȷ֪�źŵļ��
%%%%input��               T------------�۲�ʱ��
%%%%                      snr----------���ջ��˵�����ȣ�����ǰ��
%%%%                      pf-----------�޶����龯���� 
%%%%output:               pd1_est------��ͨ��ȷ֪�źŵ���ȷ������
%%%%                      pdarray_est--��ͨ��ȷ֪�źŵ���ȷ������
%%%%                      distance_est--���Ƶ�Ŀ�����
%%%%                                    ����ú���matchedfilter�����һ����Ԫ�����ź���Է����źŵ��ӳ�ʱ��Ӷ������Ŀ�����
%%%%                      doa_est-------���Ƶ�Ŀ�귽λ
%%%                       Pbeam---------�ռ��ף�����ͼ��
%%%%by tanjunhong,7.5,2011��the whole day, from moring to 20:30,the run this function for one more hour��
%%%%modified on 

% T=0.01;snr=-15;pf=0.01;  % test this function   M=10,50

%%%%%%%%%%%%%%%%��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%�źŲ���
fs=4000;                  %�Խ������ݵĲ���Ƶ��
ts=1/fs;                  %����ʱ����
N=T/ts;                   %�۲����ݸ���
t=(1:N)*ts;               %����ʱ�̵���ʸ��
c=1500;                   %����
A=1;                      %�źŷ���
fc=250;                   %�ź�Ƶ��
lamda=c/fc;               %�źŲ���
phi=pi/4;                 %�ź���λ
s=A*sin(2*pi*fc*t+phi);   %���ź�ʸ��------����ȷ֪�ź�

%%��������
snr=10^(snr/10);          %�����
Ps=A^2/2;                 %�źŵ�ƽ������
Pn=Ps/snr;                %������ƽ������
sigma=sqrt(Pn);           %�����ı�׼���ֵΪ0��

%%Ŀ�����
distance=1200;            %Ŀ�����һ����Ԫ����
doa=50;                   %Ŀ�귽λ�������з���ļнǣ�-90:90��
delay0=2*distance/c;      %��һ����Ԫ�����ӳ�tau0���յ��������ŷ����ź�
delay0_N=round(delay0/ts);
%%���в���
M=40;                     %��Ԫ��  M=10
d=lamda/2;                %���ڻ�Ԫ���
tau0=(0:M-1)'*d*sin(doa*pi/180)/c;  %������Ԫ���ӳ�ʱ�䣨��Բο���Ԫ1��
atheta0=exp(-1i*2*fc*tau0);         %����ʸ��
theta=-90:90;             %����ɨ��Ƕȷ�Χ
theta_N=length(theta);    %ɨ��Ƕ���Ŀ


%%%-���۵���ȷ������------------------------------------------------------
d1=snr*N;                 %��ͨ��
darray=snr*N*M;           %��ͨ��
u0=qfuncinv(pf);          %���龯����pf����u0  
pd1=qfunc(u0-sqrt(d1));        %��ͨ�������۵���ȷ������ֵ
pdarray=qfunc(u0-sqrt(darray));%��ͨ�������۵���ȷ������


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  ȷ����������     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G01=u0*sqrt(d1);           %��ͨ��
G0array=u0*sqrt(darray);   %��ͨ��


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   Monte Carlo  ����   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pbeam_temp=0;
pd1_est=0;
pdarray_est=0;
pf1_est=0;
pfarray_est=0;
doa_est=0;
distance_est=0;
N_exp=50;                %�����ʵ�����

for ii=1:N_exp
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % %%%%%%%%%%%H1�����£����źţ�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r=atheta0*s+sigma*randn(M,N);
    %%%% ��ͨ�����
%     G_1_1=r(1,:)*s'/Pn;    
    G_1_1=(s+sigma*randn(1,N))*s'/Pn;                        %��ͨ�����ͳ����
    if G_1_1>G01                               %�Ƚ��о�
        pd1_est=pd1_est+1;
    end
    %%%%�����źż��,�����ź�ʱ������Ŀ��������distance_est��Ŀ�귽λ����doa_est
    for k=1:theta_N
        tau=(0:M-1)'*d*sin(theta(k)*pi/180)/c; %������Ԫ���ӳ�ʱ�䣨��Բο���Ԫ1��
        atheta=exp(-1i*2*fc*tau);              %����ʸ��M*1
        beam_out(k,:)=(atheta')*r;         %�����γɣ��Ը�����Ԫ�Ľ����źŽ����ӳ����
        Pbeam(k)=beam_out(k,:)*beam_out(k,:)'; %�ռ��ף�����ͼ 
    end
    Pbeam_temp=Pbeam+Pbeam_temp;
    [Pmax,theta_index]=max(Pbeam);
    doa_est=doa_est+theta(theta_index);        %doa����
    
    G_array_1=beam_out(theta_index,:)*s'/Pn;   %���м��ͳ����
    if G_array_1>G0array                       %�Ƚ��о�
        pdarray_est=pdarray_est+1;              
    end
    
                                               %Ŀ��������
    r1=[zeros(1,delay0_N),r(1,:)];
    smf1=matchedfilter(r1,s,fs);
    [smfmax,delay_index]=max(smf1);
    distance_est=distance_est+delay_index*ts*c/2;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
% % %%%%%%%%%%%H0�����£����źţ�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r=sigma*randn(M,N);
    %%%% ��ͨ�����
    G_1_0=r(1,:)*s'/Pn;                        %��ͨ�����ͳ����
    if G_1_0>G01                               %�Ƚ��о�
        pf1_est=pf1_est+1;
    end
    %%%%�����źż��,�����ź�ʱ������Ŀ��������distance_est��Ŀ�귽λ����doa_est
    for k=1:theta_N
        tau=(0:M-1)'*d*sin(theta(k)*pi/180)/c;  %������Ԫ���ӳ�ʱ�䣨��Բο���Ԫ1��
        atheta=exp(-1i*2*fc*tau);               %����ʸ��M*1
        beam_out0(k,:)=(atheta')*r;              %�����γɣ��Ը�����Ԫ�Ľ����źŽ����ӳ����
        Pbeam0(k)=beam_out0(k,:)*beam_out0(k,:)';  %�ռ��ף�����ͼ 
    end
    Pbeam_temp0=Pbeam0+Pbeam_temp;
    [Pmax0,theta_index0]=max(Pbeam0);
    G_array_0=beam_out0(theta_index0,:)*s'/Pn;   %���м��ͳ����
    if G_array_0>G0array                       %�Ƚ��о�
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



   


