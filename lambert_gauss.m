function [x, iter] = lambert_gauss(T, lambda)

% manipulate input
tol = 1e-08;

% constants
l = (1-lambda)^2/(4*lambda);
m = T^2/(64*lambda^3);

% initial value
x = 0; %for parabola
% x = l/(1+2*l); % for elipse

iter = 0;
while true
    % evaluate continued fraction (top-down)
    odd_fn = @(n) ((n+4)*(n+7))/((2*n+5)*(2*n+7));
    even_fn = @(n) ((n-1)*(n+2))/((2*n+5)*(2*n+7));
    
    an = [];
    bn = [];
    u = [];
    del = [];
    c = [];

    % initial value
    del(end+1) = 1;
    an(end+1) = 2/35*x^2;
    bn(end+1) = 1+2/35*x;
    u(end+1) = an(end)/bn(end);
    c(end+1) = u(end);

    iter_xi = 1;
    while true
        if rem(iter_xi, 2) == 1
            an(end+1) = odd_fn(iter_xi)*x;
            bn(end+1) = 1;
        else
            an(end+1) = even_fn(iter_xi)*x;
            bn(end+1) = 1;
        end
        
        del(end+1) = 1/(1-(an(end)/(bn(end-1)*bn(end)))*del(end));
        u(end+1) = u(end)*(del(end)-1);
        c(end+1) = sum(u, 'double');

        % check convergence
        fprintf("iter: %d, del: %f, u: %f, c: %f, an: %f\n", ...
            iter_xi, del(end), u(end), c(end), an(end));
        fprintf("diff: %f\n", c(end)-c(end-1));
        if abs(c(end)-c(end-1)) < 1e-5
            xi = c(end);
            break;
        end

        iter_xi = iter_xi + 1;
        if iter_xi > 1000
            error("Inf loop occurs in calculating fraction");
        end
    end

    % compute coefficient of cubic equation
    h = m/(5/6 + l + xi);
    
    % solving cubic equaiton for y (Newton-Raphson method)
    % y_prev = 0;
    y_prev = 2/3*(1+h); % initial choice of y

    f = @(y, h) y^3 - y^2 - h*y - h/9;
    df = @(y, h) 3*y^2 - 2*y - h;
    
    iter_NR = 0;
    while true
        y = y_prev - f(y_prev, h)/df(y_prev, h);
        
        if abs(y-y_prev) < tol
            % choose positive real root
            if y < 0
                error("Not a positive real root");
            end
            break;
        end
        y_prev = y;

        iter_NR = iter_NR + 1;
        if iter_NR > 1000
            error("Inf loop occurs in solving cubic");
        end
    end

    % compute new value of x
    x_next = m/(y^2) - l;
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
