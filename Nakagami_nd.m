%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nakagami衰落函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Nakagami] = Nakagami_nd(z,Pr,m)
if (Pr == 0 || m < 0.5)
    disp('error');
else
    Nakagami = (m/Pr)^m.*z.^(m-1)/gamma(m).*exp(-m.*z/Pr);
    %Nakagami = (2*m^m*z.^(2*m-1))/(gamma(m).*Pr^m).*exp(-(m.*z)./Pr);
end