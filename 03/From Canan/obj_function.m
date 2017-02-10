function xc = obj_function(uk, K, a, b, x0, ec, t)

x = state_variables(uk, K, a, b, x0, ec, t);
x_last = x(:,end);

xc = - 1 + x_last(1) + x_last(2);
