clc; clear;
close all;

% au = 149597870.7; % 1AU (km)
% mu = 398600 / (au^3); %km^3 s^-2 % The Earth

global mu PLOT DEBUG
% mu = 4*pi^2;
PLOT = false;
DEBUG = false;

% r1 = [1, 0, 0];
% r2 = [0, 1.524, 0];
% tf = 4;
% lambmultifunc(r1, r2, tf);

% Givens
Re = 6400;
mu = 3.986e5;

chaser = OrbitElement(Re+2000, 0.002, 60, 30, 0, 0);
target = OrbitElement(Re+36000, 0.0002, 55, 35, -20, 30);

t1 = 1200:10:1400;
t2 = 1200:10:1400;
result = zeros(length(t1), length(t2));
for i = 1:length(t1)
    % calculate position vector, velocity vector
    [r1, v1, theta1] = posvel(chaser.semimajor, chaser.eccentricity, chaser.true_anomaly, t1(i));

    for j = 1:length(t2)
        [r2, v2, theta2] = posvel(target.semimajor, target.eccentricity, target.true_anomaly, t1(i)+t2(j));
        [maxN, A, E] = lambmultifunc(r1, r2, t2(j));

        J = Inf;
        for k = 1:length(A)
            if (A(k) == 0)
                continue;
            end

            vt1 = sqrt(mu/(A(k)*(1-E(k)^2))) * [-sin(theta1), E(k)+cos(theta1)];
            vt2 = sqrt(mu/(A(k)*(1-E(k)^2))) * [-sin(theta2), E(k)+cos(theta2)];
            
            dv1 = norm(vt1 - v1);
            dv2 = norm(v2 - vt2);
            
            dv = dv1 + dv2;

            if dv < J
                J = dv;
            end
        end

        if J == Inf
            J = -1;
        end
        result(i, j) = J;
        fprintf("i: %d, j: %d, J: %f\n", i, j, J);
    end
end

figure()
contour(t1,t2,result);
