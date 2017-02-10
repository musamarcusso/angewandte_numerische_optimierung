function x_pnts = state_variables(uk, K, a, b, x0, ec, t)

x_pnts = zeros(numel(x0),K);
coef = (max(t) - min(t)) / K;
for i = 1:K
    x_last = [];
    if i == 1
        x_last = x0; 
    else
        x_last = x_pnts(:,i - 1);
    end
    coll_pnts = collocation_pnts(ec(uk(i)), a, t, K, x_last);
    x_pnts(:,i) = x_last;
    for j = 1:3
        x_pnts(:,i) = x_pnts(:,i) + coef * b(j) * ec(uk(i)) * coll_pnts(:,j);
    end
end
