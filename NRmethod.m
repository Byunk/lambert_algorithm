function [x] = NRmethod(x0, f, df)
nMaxIter = 1000;
tol = 1e-08;

xPrev = x0;

iter = 0;
while true
    if nMaxIter < iter
        error('Infinite Loop in NR method');
    end

    xNext = xPrev - f(xPrev)/df(xPrev);
    if abs(xNext-xPrev) < tol
        x = xNext;
        break;
    end
    xPrev = xNext;
    iter = iter + 1;
    % fprintf("Iter : %d, x: %f\n", iter, xPrev);
end
