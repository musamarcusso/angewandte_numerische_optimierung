figure(1);
xvariables=state_variables(u_final, K, x0, ec,time_span);
time= 0:(time_span(end)/K):(time_span(end));

for l=1:K+1
    xc_last(l) = -1*obj_function(u_final, K, x0, ec, [0 time(l)]);
end
figure(2)
grid on
plot( time,[x0(:,1) xvariables(1,:)],'--r', ...
    time,[x0(:,2) xvariables(2,:)],'--b', ...
    time,xc_last,'--g');
hold on
legend('X_A','X_B','X_C')
xlabel('time')
ylabel('X_A, X_B, X_C')