function c_pnts = collocation_pnts(func, a, t, K, x_last)

coef = (max(t) - min(t)) / K;

g = sym('gk%d%d', [numel(x_last) 3]);

aux = cell(1,3);
output = [];
for i = 1:3
    aux{i} = x_last(:) - g(:,i);
    for j = 1:3
        aux{i} = aux{i} + coef * a(i,j) * func * g(:,j);
    end
    output = vertcat(output, aux{i});
end

output = solve(output);
fields = fieldnames(output);
c_pnts = zeros(numel(fields),1);
for i = 1:numel(fields)
    c_pnts(i) = subs(output.(fields{i}));
end

c_pnts = reshape(c_pnts, 3, 2)';
display(c_pnts)
