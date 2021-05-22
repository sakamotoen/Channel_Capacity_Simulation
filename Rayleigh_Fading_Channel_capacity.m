%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%瑞丽衰落下的信道容量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

gama_0 = input('请输入中止门限（dB）');
%gama_0 = 1;
u = 0:30;
len = length(u);
 %Average_SNR = 10.^((u + ((sigma.^2)*(log(10)/20)))/10);      
Average_SNR = 10.^(u/10);
%% AWGN信道容量
for i = 1:len
    C_AWGN(i) = log2(1 + Average_SNR(i));
end

figure;

plot (u,C_AWGN,'-dr');

axis([0,30,0,10]);

hold on;

xlabel('平均接收信噪比(dB)','fontsize',12);

ylabel('C/B(bit/s/Hz)','fontsize',12);

title('瑞丽衰落的容量','fontsize',14);

grid on;

%% 香农容量w TX/RX CSI

for i = 1:len
    fun = @(gama) log2(1+gama/gama_0).*Ray_nd(gama,Average_SNR(i));
    C_TXRX(i) = integral(fun,gama_0,inf);
end

plot(u,C_TXRX,'-xg');

%% 香农容量w RX CSI

for i = 1:len
    fun = @(gama) log2(1 + gama).*Ray_nd(gama,Average_SNR(i));
    C_RX(i) = integral(fun,0,inf);
end

plot(u,C_RX,'-b');

%% 最大中断容量

r_0 = 1:10;
for i = 1:10
    for j = 1:len
        fun = @(r) (1./r) .* Ray_nd(r,Average_SNR(j));
        E_r = integral(fun,r_0(i),inf);
        Ray_nd_mid = @(R) Ray_nd(R,Average_SNR(j));
        Ray_nd_r0 = integral(Ray_nd_mid,0,r_0(i));
        mid(i,j) = log2(1 + 1/(E_r))*(1-Ray_nd_r0);
        
    end
    %plot(u,mid,'+m');
end
[C_max_intrpt,index]=max(mid);
plot(u,C_max_intrpt,'-+m');

%% 零中断容量

for i = 1:len
    fun = @(r) (1./r) .* Ray_nd(r,Average_SNR(i));
    E_r_0 = integral(fun,0,inf);
    C_zreo_intrpt(i) = log2(1);
end
plot(u,C_zreo_intrpt,'*k');
legend('AWGN','TXRX CSI','RX CSI','Max intrpt','Zreo intrpt');
set(legend,'Location','northwest');