clc; clear;
close all;

% au = 149597870.7; % 1AU (km)
% mu = 398600 / (au^3); %km^3 s^-2 % The Earth

% Givens
r1 = 1.0;
r2 = 1.524;
angle = 75; %deg
Tf = 4;
MU = 4*pi^2;

% Constants
theta = deg2rad(75);
c = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(theta));
s = (r1+r2+c)/2;

%% Determination of the Minimum Transfer Time
a0 = 1.001*s/2;
[~, Lower] = wrapfunc(c, s, MU);

N = 1;
while true
    funFN = @(a) Lower.funF(a, N);
    funDFN = @(a) Lower.funDF(a, N);

    a = NRmethod(a0, funFN, funDFN);
    minTf = Lower.funTf(a, N); % Minimum transfer time for N

    if minTf > Tf
        maxN = N-1;
        fprintf("Maximum N is %d\n", maxN);
        break;
    else
        % fprintf("Minimum Transfer Time for N = %d, a = %f, t = %f\n", N, a, t);
        N = N+1;
    end
end

%% Multiple Revolution Lambert Problem Formulation (Plotting)
a = linspace(0.9, 1.5, 100000);
alpha = 2*asin(sqrt(s./(2*a)));

if (theta >= 0 && theta < pi)
    beta = 2*asin(sqrt((s-c)./(2*a)));
else
    beta = 2*pi - 2*asin(sqrt((s-c)./(2*a)));
end

funTf = @(a, alpha, beta, N) a.^(3/2).*(2*N*pi + alpha-beta ...
    - (sin(alpha)-sin(beta))) ./ sqrt(MU); % transfer time

figure();
hold on;
for N = 0:maxN
    % Upper Part
    t = funTf(a, alpha, beta, N);
    I = find(t == real(t));
    plot(a(I), t(I), 'r');

    % Lower Part
    t_minus = funTf(a, 2*pi-alpha, beta, N);
    I = find(t == real(t));
    plot(a(I), t_minus(I), 'r');
end
xlabel("Semamajor Axis, a"); ylabel("Time of Flight, tf");
ylim([0 Tf]);
title("Solutions for T >= t_minNmax");
hold off;

%% Solution to Multiple-Revolution Lambert's Problem
a0 = 1.001*s/2;
[Upper, Lower] = wrapfunc(c, s, MU);

solLambert = [];
for N = maxN:-1:1
    funGNUpper = @(a) Upper.funTf(a, N) - Tf;
    funDGNUpper = @(a) Upper.funDTf(a, N);
    aUpper = NRmethod(a0, funGNUpper, funDGNUpper);
    
    funGNLower = @(a) Lower.funTf(a, N) - Tf;
    funDGNLower = @(a) Lower.funDTf(a, N);
    aLower = NRmethod(a0, funGNLower, funDGNLower);

    if N == maxN
        if Upper.funTf(aUpper, N) >= Tf
            fprintf("upper solution at Nmax exceed given T");
            solLambert(end+1) = aUpper;
        end
        solLambert(end+1) = aLower;
    else
        solLambert(end+1) = aUpper; 
        solLambert(end+1) = aLower;
    end
    % fprintf("upper: %f, lower: %f\n", aUpper, aLower);
end

disp(solLambert);