function stop = plotfun(x, optimValues, state, K, x0, ec, t)
stop = false;

x_states = state_variables(x, K, x0, ec, t);

hold on;
grid on
subplot(1,2,1)
h1 = plot(optimValues.iteration, x_states(1), 'ro', 'linewidth', 2);
h2 = plot(optimValues.iteration, x_states(2), 'yo', 'linewidth', 2);
h3 = plot(optimValues.iteration, -optimValues.fval, 'o', 'linewidth', 2);

%legend([h1, h2, h3], 'String', {'X_A', 'X_B', 'X_C'});

hold on;
grid on
subplot(1,2,2)
plot(optimValues.iteration, norm(x,2), 'ro', 'linewidth', 2);
drawnow