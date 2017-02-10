function x_pnts = state_variables(uk, K, x0, ec, t)
% uk     Current value for the control variables given by the optimizer
% K      Number of finite elements
% x0     Initial value for the differential equation for the equality
%        constrains
% ec     Function representing all the model equations
% t      Time interval
% x_pnts Values for Xa and Xb for each of one of the K finite elements

% Coefficients of the Radau quadrature
b = [(16 - sqrt(6))/36 (16 + sqrt(6))/36 1/9]';
% Vector of the discrete points delimiting the finite elements
x_pnts = zeros(numel(x0),K);
% (T2 - T1) / K
coef = (max(t) - min(t)) / K;
% Calculating the Xak and Xbk points
for i = 1:K
    % Set the value for previous Xk point
    x_last = [];
    % For the first finite element, the previous point is the initial point
    % given, x0
    if i == 1
        x_last = x0; 
    else
        x_last = x_pnts(:,i - 1);
    end
    % Calculate the collocation points
    coll_pnts = collocation_pnts(ec(uk(i)), t, K, x_last);
    % Solve x(k) = x(k - 1) + coef * bj * collocation points
    % Remembering that Xk = [Xak Xbk]'
    x_pnts(:,i) = x_last;
    for j = 1:3
        x_pnts(:,i) = x_pnts(:,i) + coef * b(j) * ec(uk(i)) * coll_pnts(:,j);
    end
end