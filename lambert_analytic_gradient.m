function [x1, iter] = lambert_analytic_gradient(r1,...
                                                r2,...
                                                phi,...
                                                tf,...
                                                muC)

% manipulate input
tol = 1e-06;
% r1 = sqrt(r1vec*r1vec.');
% r2 = sqrt(r2vec*r2vec.');
% phi = acos(dot(r1vec, r2vec)/(r1*r2));

% calculate the bounds for the initial flight-path angle
phi = deg2rad(phi);
if (phi < pi)
    xmin = atan((cos(phi) - r1/r2)/sin(phi));
else
    xmin = -pi/2;
end
xmax = atan((sin(phi) + sqrt(2*r1/r2*(1 - cos(phi))))/(1 - cos(phi)));

% initial flight-path angle
x1 = (xmin + xmax)/2;

iter = 0;
while true
    % calculate the transfer time
    roh = tan(x1);
    vt1 = sqrt((muC*(1 - cos(phi))/(r1*(r1/r2 - cos(phi) + roh*sin(phi)))));
    vr1 = roh*vt1;
    v1 = sqrt(vt1^2 + vr1^2);
    
    a = 1/(2/r1 - v1^2/muC);        % semimajor axis
    v2 = sqrt((2/r2-1/a)*muC);
    alpha = 1/a;
%     alpha = 2/r1-v1^2/muC;
    h = r1*vt1;                      % angular momentum
    p = h^2/muC;                    % semiparameter
    e = sqrt(1-p/a);                % eccentricity
    f1 = atan2(h*vr1/muC, h*vt1/muC-1);
    f2 = f1 + phi;

    if alpha > 0 && e < 1
        % eliptic orbit
%         fprintf("Eliptic Orbit...\n");
        E1 = 2*atan(sqrt((1-e)/(1+e))*tan(f1/2)); % eccentric anomalies
        E2 = 2*atan(sqrt((1-e)/(1+e))*tan(f2/2));
        if E2 < E1
            E2 = 2*pi + E2;
        end

        n = muC^(1/2)*alpha^(3/2); % mean motion
        t = ((E2-E1) - e*(sin(E2)-sin(E1)))/n; % transfer time
    elseif alpha < 0 && e > 1
        % hyperbolic orbit
%         fprintf("Hyperbolic Orbit...\n");
        H1 = 2*atanh(sqrt((e-1)/(e+1))*tan(f1/2)); % eccentric anomalies
        H2 = 2*atanh(sqrt((e-1)/(e+1))*tan(f2/2)); 
        n = muC^(1/2)*(-alpha)^(3/2); % mean motion
        t = (e*(sinh(H2)-sinh(H1)) - (H2-H1))/n;
    elseif e == 1
        % parabolic orbit
%         fprintf("Parabolic Orbit...\n");
        n = 2*sqrt(muC/p^3);
        w1 = tan(f1/2);
        w2 = tan(f2/2);
        t = ((w2^3/3+w2) - (w1^3/3+w1))/n;
    end
    
    % check the convergence
%     fprintf("Diff: %0.10f\n", abs(t-tf));
    if abs(t-tf) < tol
        break;
    end
    
    % calculate the anlytic gradient of tf
    if alpha > 0 && e < 1
        % eliptic orbit
        G_VR1 = (r1*alpha^(1/2)*(cos(E1)-e)) / (muC^(1/2)*e);
        G_VR2 = (r2*alpha^(1/2)*(cos(E2)-e)) / (muC^(1/2)*e);
        G_alp1 = sin(E1) * ((alpha^(-1)*(cos(E1)-e)/2) + r1/e);
        G_alp2 = sin(E2) * ((alpha^(-1)*(cos(E2)-e)/2) + r2/e);
        G_alpn = 3*n*t/(2*alpha);
        
        dv1dx1N = r1*v1^3*(-sin(phi)*cos(2*x1)+(r1/r2-cos(phi))*sin(2*x1));
        dv1dx1D = 2*muC*(1-cos(phi));
        dv1dx1 = dv1dx1N/dv1dx1D; % Zeta
        dvr1dx1 = sin(x1)*dv1dx1 + v1*cos(x1);

        x2 = atan((r2/r1-1)*cot(phi/2) - r2/r1*tan(x1));
        dx2dx1 = -(r2*cos(x2)^2)/(r1*cos(x1)^2); % Gamma
        dv2dx1 = (v1/v2)*dv1dx1;
        dvr2dx1 = sin(x2)*dv2dx1 + (v2*cos(x2))*dx2dx1;

        dalpdx1 = -(2*v1)/muC*dv1dx1;

        dtdx1 = (G_VR2*dvr2dx1 - G_VR1*dvr1dx1 + (G_alp2-G_alp1-G_alpn)*dalpdx1)/n;
    elseif alpha < 0 && e > 1
        % hyperbolic orbit
        G_VR1 = (r1*(-alpha)^(1/2)*(e-cosh(H1))) / (muC^(1/2)*e);
        G_VR2 = (r2*(-alpha)^(1/2)*(e-cosh(H2))) / (muC^(1/2)*e);
        G_alp1 = sinh(H1) * ((e-cosh(H1))/(2*alpha) - r1/e);
        G_alp2 = sinh(H2) * ((e-cosh(H2))/(2*alpha) - r2/e);
        G_alpn = 3*n*t/(2*alpha);
        
        dv1dx1N = r1*v1^3*(-sin(phi)*cos(2*x1)+(r1/r2-cos(phi))*sin(2*x1));
        dv1dx1D = 2*muC*(1-cos(phi));
        dv1dx1 = dv1dx1N/dv1dx1D; % Zeta
        dvr1dx1 = sin(x1)*dv1dx1 + v1*cos(x1);

        x2 = atan((r2/r1-1)*cot(phi/2) - r2/r1*tan(x1));
        dx2dx1 = -(r2*cos(x2)^2)/(r1*cos(x1)^2); % Gamma
        dv2dx1 = (v1/v2)*dv1dx1;
        dvr2dx1 = sin(x2)*dv2dx1 + (v2*cos(x2))*dx2dx1;

        dalpdx1 = -(2*v1)/muC*dv1dx1;

        dtdx1 = (G_VR2*dvr2dx1 - G_VR1*dvr1dx1 + (G_alp2-G_alp1-G_alpn)*dalpdx1)/n;
    elseif e == 1
        q = p/(1+e);
        
        dv1dx1N = r1*v1^3*(-sin(phi)*cos(2*x1)+(r1/r2-cos(phi))*sin(2*x1));
        dv1dx1D = 2*muC*(1-cos(phi));
        dv1dx1 = dv1dx1N/dv1dx1D; % Zeta

        dhdx1 = r1*(dv1dx1*cos(x1) - v1*sin(x1)); % eta
        dedx1 = (p*v1)/muC*dv1dx1; % epsilon
        dqdx1 = -p/4*dedx1 + h/muC*dhdx1; % Theta

        df1dx1 = 2*((sin(2*x1-f1)-sin(f1))/v1*dv1dx1 + cos(2*x1-f1));
        dw1dx1 = (w1^2+1)/2*df1dx1; %Omega1

        dx2dx1 = -(r2*cos(x2)^2)/(r1*cos(x1)^2); % Gamma
        dw2dx1 = (w2^2+1)*((sin(2*x2-f2)-sin(f2))/v2*v1/v2*dv1dx1 + cos(2*x2-f2)*dx2dx1); % Omega2

        dlambdx1 = -dedx1/2;

        dtdx1 = (2*h/q*dqdx1-dhdx1)*t/h + 2*q^2/h*((1+w2^2)*dw2dx1-(1+w1^2)*dw1dx1 ...
                - 2*((w2^3/3+w2^5/5)-(w1^3/3+w1^5/5))*dlambdx1);
    end
    
    % calculate flight-path increment
    dx1 = (tf-t)/dtdx1;
%     fprintf("iter: %d, path-angle: %0.10f, delta_path_angle: %0.10f," + ...
%         "flight-time: %f, e: %f\n", iter, rad2deg(x1), rad2deg(dx1), t, e);
    
    % update flight-path angle
    % interpolation parameters at the bound
    alpha_min = 0.95;
    alpha_max = 0.85;
    if x1+dx1 >= xmin && x1+dx1 <= xmax
        x1 = x1+dx1;
    elseif x1+dx1 <= xmin
%         fprintf("Updated with alpha min\n");
        x1 = alpha_min*xmin + (1-alpha_min)*x1;
    elseif x1+dx1 >= xmax
%         fprintf("Updated with alpha max\n");
        x1 = alpha_max*xmax + (1-alpha_max)*x1;
    end

    
    iter = iter + 1;
    if iter > 100
        error("INF loop occurs");
    end
end

end

