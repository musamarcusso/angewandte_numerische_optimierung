function c_pnts = collocation_pnts(func, t, K, x_last)
% func   Model equations evaluated for the control variables given by the
%        optimizer
% t      Time interval 
% K      Number of finite elements
% x_last Value of X for the previous finite element
% c_pnts Collocation point values

% (T2 - T1) / K
coef = (max(t) - min(t)) / K;

% Creating the symbolic variables for the collocation points to be
% calculated, one line for each discretized variable Xa and Xb
% Collocation points for Xa => [ gk11, gk12, gk13]
% Collocation points for Xb => [ gk21, gk22, gk23]
g = sym('gk%d%d', [numel(x_last) 3]);

% Coefficients for Radau quadrature
a = [ (88 - 7*sqrt(6))/360  (296 - 169*sqrt(6))/1800  (-2 + 3*sqrt(6))/255; ...
    (296 + 169*sqrt(6))/1800  (88 + 7*sqrt(6))/360    (-2 - 3*sqrt(6))/255; ...
    (16 - sqrt(6))/36     (16 + sqrt(6))/36   1/9];

% Cell structure for the equations regarding all the collocation points for 
% all state variables in the current iteration
aux = cell(1,3);

% Matrix to store the completed equations
output = [];

% Calculation of gk (Eq. 3a-c)
% Substracting gk on both sides of the equation to adequate the argument
% for the solving method 
for i = 1:3
    aux{i} = x_last(:) - g(:,i);
    for j = 1:3
        aux{i} = aux{i} + coef * a(i,j) * func * g(:,j);
    end
    output = vertcat(output, aux{i});
end

% Solve the system of linear equations
output = solve(output);
fields = fieldnames(output);

c_pnts = zeros(numel(fields),1);
for i = 1:numel(fields)
    % Convert the solutions given by solve from symbolic to numeric
    % structure
    c_pnts(i) = subs(output.(fields{i}));
end

% Reshape the output matrix to the form 
% gka = [ gk11, gk12, gk13]
% gkb = [ gk21, gk22, gk23]
c_pnts = reshape(c_pnts, 3, 2)';
