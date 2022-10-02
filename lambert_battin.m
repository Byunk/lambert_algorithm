function [x, iter] = lambert_battin(T, lambda)

% manipulate input
tol = 1e-08;

% Constants
l = ((1-lambda)/(1+lambda))^2;
m = T^2 / (1+lambda)^6;

% initial value
x = 0; %for parabola
% x = l; %for ellipse

x_prev = 0;
iter = 0;
while true
    % evaluate continued fraction (improvement in the convergence)
    aNn = [];   aDn = [];
    an = [];    bn = [];
    u = [];
    del = [];
    c = [];

    % initial value
    eta = (sqrt(1+x)-1)/(sqrt(1+x)+1);

    del(end+1) = 1;
    aNn(end+1) = 1;
    aDn(end+1) = 5;
    an(end+1) = aNn(end)/aDn(end);
    bn(end+1) = 1;
    u(end+1) = an(end)/bn(end);
    c(end+1) = u(end);

    xi_prev = 0;
    iter_xi = 2;
    while true
        if iter_xi == 2
            aNn(end+1) = 9;
            aDn(end+1) = 35;
            bn(end+1) = 1;
        else
            aNn(end+1) = aNn(end)+2*iter_xi+1;
            aDn(end+1) = aDn(end)+4*(2*iter_xi+1);
            bn(end+1) = 1;
        end
        
        an(end+1) = -aNn(end)/aDn(end)*eta;
        del(end+1) = 1/(1-(an(end)/(bn(end-1)*bn(end)))*del(end));
        u(end+1) = u(end)*(del(end)-1);
        c(end+1) = sum(u, 'double');

        xi = (1/(8*(sqrt(1+x)+1)))*(3+c(end)/(1+eta*c(end)));

        % check convergence
%         fprintf("iter: %d, del: %f, u: %f, c: %f, an: %f\n", ...
%             iter_xi, del(end), u(end), c(end), an(end));
        if abs(xi-xi_prev) < tol
            break;
        end

        xi_prev = xi;
        iter_xi = iter_xi + 1;
        if iter_xi > 10000
            error("Inf loop occurs in calculating fraction");
        end
    end

    % compute coefficient of cubic equation
    hN1 = (l+x)^2 * (1+xi*(1+3*x));
    hD1 = (1+2*x+l) * (3+x * (4 * xi+1));
    h1 = hN1/hD1;

    hN2 = m * (1+xi*(x-l));
    hD2 = (1+2*x+l) * (3+x * (4 * xi+1));
    h2 = hN2/hD2;
    
    % solving cubic equaiton for y (Continued Fraction)
    B = (27*h2)/(4*(1+h1)^3);
    u = -B/(2*(sqrt(1+B)+1));

    K_I = []; % compute K(u)

    odd_fn = @(n) (2*(3*n+2)*(6*n+1)) / (9*(4*n+1)*(4*n+3));
    even_fn = @(n) (2*(3*n+1)*(6*n-1)) / (9*(4*n-1)*(4*n+1));
    
    y_prev = 0;
    iter_cubic = 1;
    while true
        if rem(iter_cubic, 2) == 1
            K_I(end+1) = odd_fn(iter_cubic)*u;
        else
            K_I(end+1) = even_fn(iter_cubic)*u;
        end

        recur_sum = 1;
        for i = length(K_I):-1:1
            recur_sum = 1 - K_I(i) / recur_sum;
        end
        
        K = (1/3) / recur_sum;
        y = ((1+h1)/3) * (2+ sqrt(1+B) / (1-2*u*K^2));

        if abs(y-y_prev) < tol
            break;
        end
        y_prev = y;

        iter_cubic = iter_cubic + 1;
        if iter_cubic > 1000
            error("Inf loop occurs in solving cubic");
        end
    end

    y = ((1+h1)/3) * (2+ sqrt(1+B) / (1-2*u*K^2));

    % compute new value of x
    x_next = 1/2*(sqrt((1-l)^2 + 4*m/y^2) - (1+l));
%     fprintf("iter: %d, x: %f\n", iter, x_next);

    % check convergence
    if abs(x-x_next) < tol
        fprintf("x converged with %d iterations\n", iter);
        x = x_next;
        break;
    end
    x = x_next;
    
    iter = iter+1;
    if iter > 1000
        error("Inf loop occurs");
    end
end

end