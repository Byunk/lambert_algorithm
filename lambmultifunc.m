function [maxN, solLambert, eLambert] = lambmultifunc(r1, r2, Tf)

% Solve multi revolution Lambert's problem

% Prussing's method

% input

% r1     = initial chaser position vector (km)
% r2     = initial target position vector (km)
% tf     = transfer time (minutes)

% output

% maxN   = maximum number of revolutions
% solLambert  = semimajor axes (2N+1)
% eLamber     = eccentricities (2N+1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global mu PLOT DEBUG

% au = 149597870.7; % 1AU (km)
% r1 = r1./au;
% r2 = r2./au;
% Tf = Tf/(60*24*365); % minutes to year

% true anomaly
theta = acos(dot(r1, r2)/(norm(r1)*norm(r2)));

% chord
c = norm(r2-r1);

% semi perimeter
s = (norm(r1)+norm(r2)+norm(c))/2;

% predefined functions class
[Upper, Lower] = wrapfunc(c, s, mu);

%% Determination of the Minimum Transfer Time
a0 = 1.001*s/2;
N = 1;
while true
    funFN = @(a) Lower.funF(a, N);
    funDFN = @(a) Lower.funDF(a, N);
    a = NRmethod(a0, funFN, funDFN);

    minTf = Lower.funTf(a, N); % Minimum transfer time for N
    
    if minTf > Tf
        maxN = N-1;
        if DEBUG
            fprintf("minTf: %f, Tf: %f\n", minTf, Tf);
            fprintf("Maximum N is %d\n", maxN);
        end
        break;
    else
        if DEBUG
            fprintf("Minimum Transfer Time for N = %d, a = %f, minTf = %f, " + ...
                "Tf = %f\n", N, a, minTf, Tf);
        end
        N = N+1;
    end
end

%% Multiple Revolution Lambert Problem Formulation (Plotting)
% a = 0.9:0.0001:1.5;
a = 0: 0.001:1.5;

alpha = 2*asin(sqrt(s./(2*a)));
if (theta >= 0)
    beta = 2*asin(sqrt((s-c)./(2*a)));
else
    beta = 2*pi - 2*asin(sqrt((s-c)./(2*a)));
end

% transfer time
funTf_lower = @(N) a.^(3/2).*(2*N*pi + alpha-beta - (sin(alpha)-sin(beta))) ./ sqrt(mu);
funTf_upper = @(N) a.^(3/2).*(2*N*pi + 2*pi-alpha-beta - (sin(2*pi-alpha)-sin(beta))) ./ sqrt(mu);

if PLOT
    figure();
    hold on;
    for N = 0:maxN
        % Lower part
        t = funTf_lower(N);
        I = find(t == real(t));
        plot(a(I), t(I), 'r');
    
        % Upper Part
        t = funTf_upper(N);
        I = find(t == real(t));
        plot(a(I), t(I), 'r');
    end

    xlabel("Semamajor Axis, a"); ylabel("Time of Flight, tf");
    title("Solutions for T >= t_minNmax");
    hold off;
end

%% Solution to Multiple-Revolution Lambert's Problem
a0 = 1.001*s/2;
funEccentricity = @(alpha, beta) sqrt(1 - 4*(s-norm(r1))*(s-norm(r2))/c^2*sin((alpha+beta)/2)^2);

solLambert = zeros(1, 2*N+1);
eLambert = zeros(1, 2*N+1);
for N = 0:maxN
    if N == 0
        funGNUpper = @(a) Upper.funTf(a, N) - Tf;
        funDGNUpper = @(a) Upper.funDTf(a, N);
        try
            aUpper = NRmethod(a0, funGNUpper, funDGNUpper);
        catch
            continue;
        end

        solLambert(2*N+1) = aUpper;
        eLambert(2*N+1) = funEccentricity(Upper.alpha(aUpper), Upper.beta(aUpper));
    elseif N == maxN
        funGNLower = @(a) Lower.funTf(a, N) - Tf;
        funDGNLower = @(a) Lower.funDTf(a, N);
        aLower = NRmethod(a0, funGNLower, funDGNLower);
    
        try
            funGNUpper = @(a) Upper.funTf(a, N) - Tf;
            funDGNUpper = @(a) Upper.funDTf(a, N);
            aUpper = NRmethod(a0, funGNUpper, funDGNUpper);
        catch
            if DEBUG
                fprintf("aUpper not converged\n");
            end
            aUpper = -1;
        end

        if aUpper == -1 || Upper.funTf(aUpper, N) >= Tf
            if DEBUG
                fprintf("upper solution at Nmax exceed given T\n");
            end
            solLambert(2*N) = aLower;
            eLambert(2*N) = funEccentricity(Lower.alpha(aLower), Lower.beta(aLower));

            solLambert = solLambert(1:end-1);
            eLambert = eLambert(1:end-1);
        else
            solLambert(2*N) = aLower;
            eLambert(2*N) = funEccentricity(Lower.alpha(aLower), Lower.beta(aLower));
            solLambert(2*N+1) = aUpper;
            eLambert(2*N+1) = funEccentricity(Upper.alpha(aUpper), Upper.beta(aUpper));
        end
    else
        funGNLower = @(a) Lower.funTf(a, N) - Tf;
        funDGNLower = @(a) Lower.funDTf(a, N);
        aLower = NRmethod(a0, funGNLower, funDGNLower);

        solLambert(2*N) = aLower;
        eLambert(2*N) = funEccentricity(Lower.alpha(aLower), Lower.beta(aLower));
    
        funGNUpper = @(a) Upper.funTf(a, N) - Tf;
        funDGNUpper = @(a) Upper.funDTf(a, N);
        aUpper = NRmethod(a0, funGNUpper, funDGNUpper);

        solLambert(2*N+1) = aUpper;
        eLambert(2*N+1) = funEccentricity(Upper.alpha(aUpper), Upper.beta(aUpper));
    end

    if DEBUG
%         fprintf("upper: %f, lower: %f\n", aUpper, aLower);
    end
end

if DEBUG
    disp("Solution (semi-major axis)"); disp(solLambert);
    disp("Solution (eccentricity)"); disp(eLambert);
end

end

