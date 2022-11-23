clc; clear;
close all;

global mu PLOT DEBUG
mu = 3.986e5 * 3600; %km^3min^-2
PLOT = false;
DEBUG = true;

% r1 = [1, 0, 0];
% r2 = [0, 1.524, 0];
% tf = 4;
% lambmultifunc(r1, r2, tf);

% Givens
% km, deg
Re = 6400;
chaser = OrbitElement(Re+2000, 0.002, 60, 30, 0, 0);
target = OrbitElement(Re+36000, 0.0002, 55, 35, -20, 30);

% min
t1 = 0:1:1400;
t2 = 0:1:1400;
result = zeros(length(t1), length(t2));
for i = 1:length(t1)
    % calculate position vector, velocity vector
    % [r1, v1, theta1] = posvel(chaser.semimajor, chaser.eccentricity, chaser.true_anomaly, t1(i));
    r1 = chaser.pos(t1(i));
    v1 = chaser.vel(t1(i));
    theta1 = chaser.trueanomalyAfter(t1(i));

    for j = 1:length(t2)
        % [r2, v2, theta2] = posvel(target.semimajor, target.eccentricity, target.true_anomaly, t1(i)+t2(j));
        r2 = target.pos(t2(i));
        v2 = target.vel(t2(i));
        theta2 = target.trueanomalyAfter(t2(i));

        [maxN, A, E] = lambmultifunc(r1, r2, t2(j));

        J = Inf;
        for k = 1:length(A)
            if (A(k) == 0) || A(k) ~= real(A(k)) || E(k) ~= real(E(k))
                continue;
            end

            vt1 = sqrt(mu/(A(k)*(1-E(k)^2))) * [-sin(theta1), E(k)+cos(theta1)];
            vt2 = sqrt(mu/(A(k)*(1-E(k)^2))) * [-sin(theta2), E(k)+cos(theta2)];
            fprintf("vt1: %p, v1: %p\n", norm(vt1), norm(v1));
            
            dv1 = norm(vt1 - v1);
            dv2 = norm(v2 - vt2);
            
            dv = dv1 + dv2;

            if dv < J
                J = dv;
            end
        end

        if J == Inf
            J = 0;
        end
        result(i, j) = J;
        fprintf("i: %d, j: %d, J: %f\n", i, j, J);
    end
    
end

save('result.mat','result');

%%
figure()
contour(t2, t1, result)
