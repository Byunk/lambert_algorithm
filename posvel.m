function [r, v, theta1] = posvel(a, e, theta0, t)

% Calculate position vector after given time given orbit element

% input
% a         = semimajor axis
% e         = eccentricity
% theta0    = initial true anomaly
% tp        = parking time (seconds)

% output

% r         = position vector after tp seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global mu

% a = elem.semimajor;
% e = elem.eccentricity;
% theta0 = elem.true_anomaly;

theta0 = deg2rad(theta0);

E0 = 2*atan(tan(theta0/2) * sqrt((1-e)/(1+e))); % Eccentric Anomaly
M0 = (E0 - e*sin(E0)); % Mean Anomaly
M1 = M0 + sqrt(mu/a^2) * t;

% function define for Newton-Raphson method
f = @(E) M1 - (E - e*sin(E));
df = @(E) -(1 - e*cos(E));
E1 = NRmethod(E0, f, df);

theta1 = 2*atan(tan(E1/2) * sqrt((1+e)/(1-e))); % true anoamly after tp
r = a*(1-e^2)/(1+e*cos(theta1)) * [cos(theta1), sin(theta1)];
v = sqrt(mu/(a*(1-e^2))) * [-sin(theta1), e+cos(theta1)];
end

