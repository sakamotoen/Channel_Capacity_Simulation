clc;close all;clear
Nsample = 100000;
m = 2;%Nakagami-m,m=2 
Average_SNR_dB = 0:1:30;
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
pdfx = @(x,mu,m) (m/mu)^m*x.^(m-1)/gamma(m).*exp(-m*x/mu);
for index_snr0 = 1:length(Average_SNR)
    mu = Average_SNR(index_snr0);
    gamma0 = 1e-3;
    abserr = 1;
    while(abserr>1e-2)
        gamma0 = gamma0 + 1e-3;
        abserr = abs(1 - integral(@(x)(1/gamma0 - 1./x).*pdfx(x,mu,2),gamma0,Inf,'AbsTol',1e-3,'RelTol',1e-2)); % integral Introduced in R2012a
    end
    CSITR_gamma0(index_snr0) = gamma0;
end
% search TCI_gamma0
for index_snr1 = 1:length(Average_SNR)
    Range_gamma0 = 0:1:200;
    accuracy = 1;
    mu = Average_SNR(index_snr1);
    while(accuracy >10^-2)
        accuracy = accuracy/10;
        E_gamma0 = zeros(size(Range_gamma0));
        P1 = zeros(size(Range_gamma0));% P(gamma >= gamma0 )
        C_Pout = zeros(size(Range_gamma0));
        for index_gamma = 1:length(Range_gamma0)
            gamma0 = Range_gamma0(index_gamma);
            E_gamma0 = integral(@(x)(1./x).*pdfx(x,mu,2),gamma0,Inf,'AbsTol',1e-3,'RelTol',1e-2);%E_gamma0{1/gamma} (4.20)
            Pout = integral(@(x)pdfx(x,mu,2),0,gamma0,'AbsTol',1e-3,'RelTol',1e-2);            
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
E1 = zeros(size(Average_SNR));%E{1/gamma}(4.18)
for i = 1:Nsample
    gamma_i = gamrnd(m,Average_SNR/m);
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
plot(Average_SNR_dB,AWGN_Capacity,'k-*')
plot(Average_SNR_dB,Shannon_Capacity_CSITR,'b-')
plot(Average_SNR_dB,Shannon_Capacity_CSIR,'r--')
plot(Average_SNR_dB,Maximum_Outage_Capcity,'g-+')
plot(Average_SNR_dB,Zero_Outage_Capcity,'r-o')
grid on
legend1 = legend('AWGN Channel Capacity','Shannon Capacity w TX/RX CSI(4.9)','Shannon Capacity w RX CSI(4.5)','Maximum Outage Capacity (4.22)','Zero-Outage Capacity (4.18)');
set(legend1,'Location','northwest');
%axis([0,30 0,10])
xlabel('Average SNR(dB)')
ylabel('C/B (Bits/Sec/Hz)')
title('Figure 4.8:Capacity in Nakagami Fading(m=2)')