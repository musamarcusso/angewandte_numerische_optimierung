function xc = obj_function(uk, K, x0, ec, t)
% uk  Input vector for the control variables
% K   Number of finite elements
% x0  Initial values of the model
% ec  Function representing all the model equations
% t   Time interval

% Calculate the values for [Xak Xbk]'
x = state_variables(uk, K, x0, ec, t);
% Values of Xa and Xb for the last finite element
x_last = x(:,end);
% Calculation of Xc
xc = - 1 + x_last(1) + x_last(2);
