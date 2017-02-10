%% Programming exercise #3
%  Dynamic optimization

clc
clear all
close all

time_span = [0 1];
K = 5;

xa0 = 1;
xb0 = 0;

ec = @(u)[-u 10*u; u, -9*u-1];
x0 = [xa0 xb0];

% Coefficients for Radau quadrature
a = [ (88 - 7*sqrt(6))/360  (296 - 169*sqrt(6))/1800  (-2 + 3*sqrt(6))/255; ...
    (296 + 169*sqrt(6))/1800  (88 + 7*sqrt(6))/360    (-2 - 3*sqrt(6))/255; ...
    (16 - sqrt(6))/36     (16 + sqrt(6))/36   1/9];

b = [(16 - sqrt(6))/36 (16 + sqrt(6))/36 1/9]';
for k = K
    
    u = struct();
    u.lb = zeros(K,1);
    u.ub = ones(K,1);
    u.u0 = 0.5 * ones(K,1);
    
    opt_options = optimset('Algorithm', 'interior-point', 'Display', 'iter', ...
        'PlotFcns',{@optimplotx,@optimplotfval});
    
%     opt_options = optimset('Algorithm', 'interior-point', 'Display', 'iter', ...
%         'OutputFcn', @(x, optimValues, state)plotfun(x, optimValues, state, K, a, b, x0, ec, time_span));
%     
    
    [u_final, fval, exitflag] = fmincon(@(uk)obj_function(uk, k, a, b, x0, ec, time_span), ...
        u.u0, [], [], [], [], u.lb, u.ub, [], opt_options);
    
%     figure
%     subplot(1,2,2)
%     t = linspace(min(time_span), max(time_span), k);
%     bar(t,u_final, 'histc');
    
end
run 
run
