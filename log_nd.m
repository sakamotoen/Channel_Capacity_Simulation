%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%该函数为对数正态分布概率密度函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[p] = log_nd(x,mu,sigma)
zeta = 10/log(10);
if x == 0
    disp('error');
else
    p = zeta./(sqrt(2*pi)*sigma.*x).*exp(-(10*log10(x)-mu).^2./(2*sigma^2));
    
end