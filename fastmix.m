function [xt] = fastmix(x0,K,rho)
% FastMix algorithm for average consensus
x_1 = x0;
eta = (1-sqrt(1-rho^2))/(1+sqrt(1-rho^2));
global W;
for k=1:K
    xt = (1+eta)*(W*x0) - eta*x_1;
    x_1 = x0;
    x0 = xt;
end
end

