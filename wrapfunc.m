function [Upper, Lower] = wrapfunc(c, s, MU)
%% Lower Portion
alpha = @(a) 2*asin(sqrt(s/(2*a)));
beta = @(a) 2*asin(sqrt((s-c)/(2*a)));
zeta = @(a) alpha(a)-beta(a);
eta = @(a) sin(alpha(a))-sin(beta(a));

Lower.funF = @(a, N) (6*N*pi + 3*zeta(a) - eta(a)) ...
                * (sin(zeta(a)) + eta(a)) - 8*(1-cos(zeta(a)));
Lower.funDF = @(a, N) ((6*N*pi + 3*zeta(a) - eta(a))*(cos(zeta(a)) + cos(alpha(a)))...
            + (3-cos(alpha(a)))*(sin(zeta(a))+eta(a))-8*sin(zeta(a))) * (-1/a * tan(alpha(a)/2))...
            + ((6*N*pi + 3*zeta(a) - eta(a)) * (-cos(zeta(a))-cos(alpha(a)))...
            + (-3-cos(beta(a))) * (sin(zeta(a))+eta(a)) + 8*sin(zeta(a))) * (-1/a * tan(beta(a)/2));
Lower.funTf = @(a, N) a^(3/2) * (2*N*pi + zeta(a) - eta(a)) / sqrt(MU);
Lower.funDTf = @(a, N) 1/2*(a/MU)^(1/2) / (sin(zeta(a)) + eta(a)) * Lower.funF(a, N);

%% Upper Portion
alpha = @(a) 2*pi - 2*asin(sqrt(s/(2*a)));
beta = @(a) 2*asin(sqrt((s-c)/(2*a)));
zeta = @(a) alpha(a)-beta(a);
eta = @(a) sin(alpha(a))-sin(beta(a));

Upper.funF = @(a, N) (6*N*pi + 3*zeta(a) - eta(a)) ...
                * (sin(zeta(a)) + eta(a)) - 8*(1-cos(zeta(a)));
Upper.funDF = @(a, N) ((6*N*pi + 3*zeta(a) - eta(a))*(cos(zeta(a)) + cos(alpha(a)))...
            + (3-cos(alpha(a)))*(sin(zeta(a))+eta(a))-8*sin(zeta(a))) * (-1/a * tan(alpha(a)/2))...
            + ((6*N*pi + 3*zeta(a) - eta(a)) * (-cos(zeta(a))-cos(alpha(a)))...
            + (-3-cos(beta(a))) * (sin(zeta(a))+eta(a)) + 8*sin(zeta(a))) * (-1/a * tan(beta(a)/2));
Upper.funTf = @(a, N) a^(3/2) * (2*N*pi + zeta(a) - eta(a)) / sqrt(MU);
Upper.funDTf = @(a, N) 1/2*(a/MU)^(1/2) / (sin(zeta(a)) + eta(a)) * Upper.funF(a, N);
end

