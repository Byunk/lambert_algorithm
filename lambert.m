function [x, iter] = lambert(T, lambda, method)

switch method
    case 'Gauss'
        [x, iter] = lambert_gauss(T, lambda);
    case 'Battin'
        [x, iter] = lambert_battin(T, lambda);
    case 'Analytic_Gradient'
        [x, iter] = lambert_analytic_gradient();
    otherwise
        error("Wrong Method");
end

end