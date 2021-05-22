clc;close all;clear
Nsample = 100000;
Average_dB_SNR = 0:2.5:30;%average dB SNR: miu_dB
sigma_dB = 8;
Average_SNR_dB = Average_dB_SNR + sigma_dB.^2*log(10)/20; %average SNR in dB:10log_10(miu) = miu_dB + sigma_dB^2ln(10)/20   (2.46)
Average_SNR = 10.^(Average_SNR_dB/10);
%(1) AWGN_Capacity
AWGN_Capacity = log2(1+Average_SNR);
% flat fading Capacity
Shannon_Capacity_CSITR = zeros(size(Average_SNR_dB));
Shannon_Capacity_CSIR = zeros(size(Average_SNR_dB));
Maximum_Outage_Capcity = zeros(size(Average_SNR_dB));
Zero_Outage_Capcity = zeros(size(Average_SNR_dB));%for Rayleigh fading: Zero_Outage_Capcity = 0;

CSITR_gamma0 = zeros(size(Average_SNR_dB));%Shannon_Capacity_CSITR: waterfilling cutoff fade depth
TCI_gamma0 = zeros(size(Average_SNR_dB));%Maximum_Outage_Capcity: truncated channel inversion cutoff fade depth

%calculate CSITR_gamma0
pdfx = @(x,mu_dB,sigma_dB) 10/log(10)/sqrt(2*pi)./sigma_dB./x.*exp(-(10*log10(x)-mu_dB).^2/2./sigma_dB^2);
for index_snr0 = 1:length(Average_dB_SNR)
    mu_dB = Average_dB_SNR(index_snr0);
    gamma0 = 1e-3;
    abserr = 1;
    while(abserr>1e-2)
        gamma0 = gamma0 + 1e-3;
        abserr = abs(1 - integral(@(x)(1/gamma0 - 1./x).*pdfx(x,mu_dB,sigma_dB),gamma0,Inf,'AbsTol',1e-3,'RelTol',1e-2)); % integral Introduced in R2012a        
    end
    CSITR_gamma0(index_snr0) = gamma0;
end
% search TCI_gamma0
for index_snr1 = 1:length(Average_SNR)
    Range_gamma0 = 1e-2:1:200;
    accuracy = 1;
    mu_dB = Average_dB_SNR(index_snr1);
    while(accuracy >10^-2)
        accuracy = accuracy/10;
        E_gamma0 = zeros(size(Range_gamma0));
        P1 = zeros(size(Range_gamma0));% P(gamma >= gamma0 )
        C_Pout = zeros(size(Range_gamma0));
        for index_gamma = 1:length(Range_gamma0)
            gamma0 = Range_gamma0(index_gamma);
            E_gamma0 = integral(@(x)1./x.*pdfx(x,mu_dB,sigma_dB),gamma0,Inf,'AbsTol',1e-3,'RelTol',1e-2);
            Pout = integral(@(x)pdfx(x,mu_dB,sigma_dB),0,gamma0,'AbsTol',1e-3,'RelTol',1e-2);
            C_Pout(index_gamma) = log2(1+1/E_gamma0)*(1-Pout);
        end
        [C_Pout_max , index_gamma] = max(C_Pout);
        snr1_gamma0 = Range_gamma0(index_gamma);
        Range_gamma0 = max(Range_gamma0(index_gamma) - accuracy,0): accuracy/10 :Range_gamma0(index_gamma) + accuracy;      
    end
    TCI_gamma0(index_snr1) = snr1_gamma0;
    Maximum_Outage_Capcity(index_snr1) = C_Pout_max;
end
%simulation:  Shannon_Capacity_CSIR,Shannon_Capacity_CSITR,Zero_Outage_Capcity
E1 = zeros(size(Average_dB_SNR));
for i = 1:Nsample
    gamma_i_dB = normrnd(Average_dB_SNR,sigma_dB);
    gamma_i = 10.^(gamma_i_dB/10);
    Shannon_Capacity_CSIR = Shannon_Capacity_CSIR + log2(1+gamma_i);
    S_gamma = zeros(size(Average_SNR_dB));
    index = find(gamma_i >= CSITR_gamma0);% P(gamma)/Average_P   (4.24)
    S_gamma(index) = 1./CSITR_gamma0(index) - 1./gamma_i(index);
    Shannon_Capacity_CSITR = Shannon_Capacity_CSITR + log2(1+ gamma_i .* S_gamma);
    E1 = E1 + 1./gamma_i;
end
Shannon_Capacity_CSIR = Shannon_Capacity_CSIR/Nsample;
Shannon_Capacity_CSITR = Shannon_Capacity_CSITR/Nsample;
E1 =  E1/Nsample;
Zero_Outage_Capcity = log2(1+1./E1);
figure
hold on
plot(Average_dB_SNR,AWGN_Capacity,'k-*')
plot(Average_dB_SNR,Shannon_Capacity_CSITR,'b-')
plot(Average_dB_SNR,Shannon_Capacity_CSIR,'r--')
plot(Average_dB_SNR,Maximum_Outage_Capcity,'g-+')
plot(Average_dB_SNR,Zero_Outage_Capcity,'r-o')
grid on
legend1 = legend('AWGN Channel Capacity','Shannon Capacity w TX/RX CSI(4.9)','Shannon Capacity w RX CSI(4.5)','Maximum Outage Capacity (4.22)','Zero-Outage Capacity (4.18)');
set(legend1,'Location','northwest');
%axis([5,30 0,14])
xlabel('Average dB SNR(dB)')
ylabel('C/B (Bits/Sec/Hz)')
title('Figure 4.6:Capacity in Log-Normal Shadowing')
