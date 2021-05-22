%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%对数正态衰落下的信道容量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

gama_0 = input('请输入中止门限（dB）');
%gama_0 = 1;
sigma = 8;
u = 5:30;                       %分贝信噪比的平均值u
len = length(u);
u_e = 10.^((u + ((sigma.^2)*(log(10)/20)))/10);           %平均信噪比gama

%% AWGN信道容量
for i = 1:len
    C_AWGN(i) = log2(1 + u_e(i));
end

figure;

plot (u,C_AWGN,'-dr');

axis([5,30,0,14]);

hold on;

xlabel('平均接收信噪比(dB)','fontsize',12);

ylabel('C/B(bit/s/Hz)','fontsize',12);

title('对数正态衰落下的信道容量','fontsize',14);

grid on;

%% 香农容量w TX/RX CSI

for i = 1:len
    fun = @(gama) log2(gama/gama_0).*log_nd(gama,u(i),sigma);
    C_TXRX(i) = integral(fun,gama_0,inf);
end

plot(u,C_TXRX,'-xg');

%% 香农容量w RX CSI

for i = 1:len
    fun = @(gama) log2(1 + gama).*log_nd(gama,u(i),sigma);
    C_RX(i) = integral(fun,0,inf);
end

plot(u,C_RX,'-b');

%% 最大中断容量

r_0 = 1:10;
for i = 1:10
    for j = 1:len
        fun = @(r) (1./r) .* log_nd(r,u(j),sigma);
        E_r = integral(fun,r_0(i),inf);
        log_nd_mid = @(R) log_nd(R,u(j),sigma);
        log_nd_r0 = integral(log_nd_mid,r_0(i),inf);
        mid(i,j) = log2(1 + 1/(E_r))*log_nd_r0;
    end
end
[C_max_intrpt,index]=max(mid);
plot(u,C_max_intrpt,'-+m');

%% 零中断容量

for i = 1:len
    fun = @(r) (1./r) .* log_nd(r,u(i),sigma);
    E_r_0 = integral(fun,0,inf);
    C_zreo_intrpt(i) = log2(1 + 1/E_r_0);
end
plot(u,C_zreo_intrpt,'-*k');
legend('AWGN','TXRX CSI','RX CSI','Max intrpt','Zreo intrpt');
