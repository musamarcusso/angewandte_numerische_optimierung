function [vub, vlb] = calc_bounds_extra(ub,lb,i,j)

f = {@(x,y)(lb(i).*y + x.*lb(j) - lb(i).*lb(j)); ...
    @(x,y)(ub(i).*y + x.*ub(j) - ub(i).*ub(j)); ...
    @(x,y)(ub(i).*y + x.*lb(j) - ub(i).*lb(j)); ...
    @(x,y)(lb(i).*y + x.*ub(j) - lb(i).*ub(j))};

comb = [lb(i) lb(j);...
    lb(i) ub(j); ...
    ub(i) lb(j); ...
    ub(i) ub(j)];
auxub = [];
auxlb = [];

for m = 1:size(comb,1)
    auxlb = [auxlb f{1}(comb(m,1), comb(m,2))];
    auxlb = [auxlb f{2}(comb(m,1), comb(m,2))];
    auxub = [auxub f{3}(comb(m,1), comb(m,2))];
    auxub = [auxub f{4}(comb(m,1), comb(m,2))];
end
vub = max(auxub);
vlb = min(auxlb);
vub = ub(i)*ub(j);
vlb = lb(i)*lb(j);