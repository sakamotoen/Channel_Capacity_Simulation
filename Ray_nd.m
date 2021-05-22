%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%瑞利衰落函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ray] = Ray_nd(z,Pr)
if Pr == 0
    disp('error');
else
    %ray = (2*z./Pr).*exp(-z.^2./Pr);
    ray = (1/Pr).*exp(-z./Pr);
end