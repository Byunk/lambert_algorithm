clear variables; close all; clc;

%% Input

%'Gauss', 'Battin', 'Analytic_Gradient'
% method = 'Gauss';

%% Test
lambda = [-0.9, -0.7, -0.5, -0.3, -0.1, 0.0, ...
    0.1, 0.3, 0.5, 0.7, 0.9];
T = [0.3, 0.5, 0.7, 0.9, 1.0, 3.0, 5.0, 7.0, 9.0, 11.0];

Re = 6.3781*10^3;
muC = 3.986*10^5;
r1 = [Re, Re, Re, 10*Re, Re];
r2 = [Re, Re/2, 2*Re, 10*Re, 10*Re];
phi = [5, 30, 55, 80, 105, 130, 155, 180, 205, 230, 255, 280, 305, 330, 355];
tf = [50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200];

%% Gauss method
result = zeros(length(lambda), length(T));
for i = 1:length(T)
    column = [];
    for j = 1:length(lambda)
        fprintf("Gauss method for T = %f, lambda = %f\n", ...
            T(i), lambda(j));
        try
            [~, iter] = lambert_gauss(T(i), lambda(j));
            result(j, i) = iter;
        catch
            result(j, i) = "err";
        end
    end
end

result_table = array2table(result, "VariableNames", ...
    string(T), "RowNames", string(lambda));
file_name = strcat("Gauss", ".csv");
writetable(result_table, file_name);
disp(result_table);

%% Battin method (Free Parameter)
result = zeros(length(lambda), length(T));
for i = 1:length(T)
    column = [];
    for j = 1:length(lambda)
        fprintf("Battin method for T = %f, lambda = %f\n", ...
            T(i), lambda(j));
        try
            [~, iter] = lambert_battin(T(i), lambda(j));
            result(j, i) = iter;
        catch
            result(j, i) = "err";
        end
    end
end

result_table = array2table(result, "VariableNames", ...
    string(T), "RowNames", string(lambda));
file_name = strcat("Battin", ".csv");
writetable(result_table, file_name);
disp(result_table);

%% Battin method (Iterations near 360 degrees)
lambda_for_iter = -0.99:0.01:-0.90;
T_for_iter = 1:2:11;

result = zeros(length(lambda_for_iter), length(T_for_iter));
for i = 1:length(T_for_iter)
    for j = 1:length(lambda_for_iter)
        fprintf("Battin method for T = %f, lambda = %f\n", ...
            T_for_iter(i), lambda_for_iter(j));
        try
            [~, iter] = lambert_battin(T_for_iter(i), ...
                lambda_for_iter(j));
            result(j, i) = iter;
        catch
            result(j, i) = "err";
        end
    end
end

result_table = array2table(result, "VariableNames", ...
    string(T_for_iter), "RowNames", string(lambda_for_iter));
disp(result_table);

%% Analytic Gradient method
test_titles = [
    "ReRe", ...
    "ReRe_2", ...
    "Re2Re", ...
    "10ReRe", ...
    "Re10Re"
    ];

for i = 1:length(r1)
    result_row = zeros(length(phi), length(tf));
    for j = 1:length(phi)
        for k = 1:length(tf)
%             fprintf("Analytic Gradient method for phi = %d, " + ...
%                 "tf = %d\n", phi(j), tf(k));
            try
                [~, iter] = lambert_analytic_gradient( ...
                    r1(i), r2(i), phi(j), tf(k), muC);
                result_row(j, k) = iter;
            catch
                result_row(j, k) = "err";
            end
        end
    end
    result_table = array2table(result_row, "VariableNames", ...
        string(tf), "RowNames", string(phi));
    disp(result_table);
    file_name = strcat("AG", test_titles(i), ".csv");
    writetable(result_table, file_name);
end

%% Comparison (Above 2 proposed methods)
phi = [5, 30, 55, 80, 105, 130, 155, 180, 205, 230, ...
    255, 280, 305, 330, 355];
tf = [50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200];

% Re, Re
result_row1 = zeros(length(phi), length(tf));
result_row2 = zeros(length(phi), length(tf));
for i=1:length(phi)
    lambda = cos(0.5*deg2rad(phi(i)));
    s = sqrt(2*Re^2*(1-cos(phi(i))));
    
    for j = 1:length(tf)
        T = sqrt(8*muC/s^3)*tf(j);
        try
            [~, iter1] = lambert_analytic_gradient( ...
                Re, Re, phi(i), tf(j), muC);
        catch
            iter1 = Inf;
        end
        
        try 
            [~, iter2] = lambert_battin(T, ...
                lambda);
        catch
            iter2 = Inf;
        end
        fprintf("iter1: %d, iter2: %d\n", iter1, iter2);
        result_row1(i, j) = iter2;
        result_row2(i, j) = iter1;
    end
end

result_table1 = array2table(result_row1, "VariableNames", ...
        string(tf), "RowNames", string(phi));
result_table2 = array2table(result_row2, "VariableNames", ...
        string(tf), "RowNames", string(phi));
disp(result_table1);
disp(result_table2);