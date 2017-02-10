%% Programming exercise #3
%  Dynamic optimization

clc
clear all
close all
% Definition of the time span where the function should be evaluated
time_span = [0 1];
% Number of finite elements for later discretisation of the equality
% constrains
K = 5;

% Initial conditions for Xa and Xb
xa0 = 1;
xb0 = 0;

% The equality constrains in matrix form
ec = @(u)[-u 10*u; u, -9*u-1];
% Initial conditions in vector form
x0 = [xa0 xb0];

% For loop built if an evaluation for different K values is needed
for k = K
    % Initial conditions for the control variables
    u = struct();
    u.lb = zeros(K,1);
    u.ub = ones(K,1);
    u.u0 = 0.5 * ones(K,1);
    % Setting the options for the optimization algorithm
    opt_options = optimset('Algorithm', 'interior-point', 'Display', 'final', ...
        'PlotFcns',{@optimplotx,@optimplotfval});
    % Applying fmincon to find the best values for u(t) for each finite
    % elements k
    [u_final, fval, exitflag] = fmincon(@(uk)obj_function(uk, k, x0, ec, time_span), ...
        u.u0, [], [], [], [], u.lb, u.ub, [], opt_options);
    
    disp(sprintf('Final value for Xc(t = 1): %f', -fval))
end

% Plot the results 
run;
