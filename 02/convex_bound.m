%% CONVEX BOUND
function [f_lb, f_ub] = convex_bound(n, m, c, H, Q, A, b, lb, ub)

f_lb = [];
f_ub = [];
% Testing the input data dimensions
if ~isequal(size(lb), [n 1]) || ~isequal(size(ub), [n 1])
    error('lb and ub should have length %d', n)
end

if ~isequal(size(Q), [m*n n])
    error('Q matrix should have size %d x %d', m*n, n)
end

if ~isequal(size(H), [n n])
    error('H matrix should have size %d x %d', n, n)
end

if ~isequal(size(c), [n 1])
    error('c vector should have size %d x %d', n, 1)
end

if ~isequal(size(b), [m 1])
    error('b vector should have size %d x %d', m, 1)
end

if ~isequal(size(A), [m n])
    error('A matrix should have size %d x %d', m, n)
end

%% Test the equality contrains' matrices symmetry
Qm = cell(m,1);
count = 1;
for i = 1:n:n*m
    Qm{count} = Q(i:(i + n - 1), :);
    if ~isequal(Qm{count}, Qm{count}')
        error('Q matrix #%d is not symmetric',i)
    end
    count = count + 1;
end

%% Compute the auxiliary variables
% Construct a binary matrix to check the positions in the Q matrices that
% have values different than 0 (to avoid repeating auxiliary variables
% creation)

w = zeros(n,n);
for i = 1:numel(Qm)
    w = w | (Qm{i} ~= 0);
end

%% Create the contrains for the QP problem
% Order the indexes of the pairs of bilinear terms by search non-zero
% elements on the matrix w
w_idx = [];
for i = 1:n
    for j = i:n % Search the elements of the upper triangle of the matrix
        if w(i,j) == 1
            w_idx = vertcat(w_idx, [i,j]);
        end
    end
end

% Save the number of auxiliary variables
nw = size(w_idx,1);

% Create the matrices for the inequality constrains of the QP problem
Aineq = [];
bineq = [];

for i = 1:nw
    % For each additional auxiliary variable there are 4 inequality
    % constrains needed
    Aaux = zeros(nw, n + nw);
    baux = zeros(nw, 1);
    
    % The inequality constrains are described by the McCormick envelopes
    % xi*lb(j) + xj*lb(i) - w_ij <= lb(i)*lb(j)
    % xi*ub(j) + xj*ub(i) - w_ij <= ub(i)*ub(j)
    % -xi*lb(j) - xj*ub(i) + w_ij <= -ub(i)*lb(j)
    % -xi*ub(j) - xj*lb(i) + w_ij <= -lb(i)*ub(j)
    
    Aaux(1, w_idx(i,1)) = lb(w_idx(i,2));
    Aaux(1, w_idx(i,2)) = lb(w_idx(i,1));
    
    Aaux(2, w_idx(i,1)) = ub(w_idx(i,2));
    Aaux(2, w_idx(i,2)) = ub(w_idx(i,1));
    
    Aaux(3, w_idx(i,1)) = -ub(w_idx(i,2));
    Aaux(3, w_idx(i,2)) = -lb(w_idx(i,1));
    
    Aaux(4, w_idx(i,1)) = -lb(w_idx(i,2));
    Aaux(4, w_idx(i,2)) = -ub(w_idx(i,1));
    
    Aaux(:, n+i) = [-1 -1 1 1]';
    
    % Add the inequality constrains for the auxiliary variable w_ij to the
    % inequality constrain matrix Aineq and vector bineq
    Aineq = vertcat(Aineq, Aaux);
    
    baux = [lb(w_idx(i,1))*lb(w_idx(i,2)); ...
        ub(w_idx(i,1))*ub(w_idx(i,2)); ...
        -ub(w_idx(i,1))*lb(w_idx(i,2)); ...
        -lb(w_idx(i,1))*ub(w_idx(i,2))];
    
    bineq = vertcat(bineq, baux);
end

% Adequate the dimension of A including the created auxiliary variables
newA = horzcat(A, zeros(size(A,1), nw));

% Modify equality constrain matrix A and vector b using the given
% coefficients from the nonlinear terms from Q matrices
for i = 1:m
    for j = 1:size(w_idx,1)
        if Qm{i}(w_idx(j,1), w_idx(j,2)) ~= 0
            if w_idx(j,1) == w_idx(j,2)
                % Coefficients of quadratic terms
                newA(i,n+j) = Qm{i}(w_idx(j,1), w_idx(j,2));
            else
                % Coefficients of bilinear terms should be multiplied by
                % factor 2
                newA(i,n+j) = 2 * Qm{i}(w_idx(j,1), w_idx(j,2));
            end
        end
    end
end

% Adequate the dimensions of H and c for the quadratic problem
newH = zeros(n+nw, n+nw);
newc = zeros(n+nw,1);

% Set the original matrix H and vector c for the variables xi
newH(1:n,1:n) = H;
newc(1:n) = c;

% Include the lower and upper bound of the auxiliary variables for the
% quadratic problem
newLB = zeros(n+nw,1);
newLB(1:n) = lb;
newUB = zeros(n+nw,1);
newUB(1:n) = ub;
% Calculate the lower and upper bounds for the auxiliary variables w_ij
for i = 1:nw
    [upper, lower] = calc_bounds(ub, lb, w_idx(i,1), w_idx(i,2));
    newLB(n+i) = lower;
    newUB(n+i) = upper;
end

%% Find the lower bound using the QP solver applying convex relaxation on original problem
opt = optimset('Display','off','Algorithm', 'active-set');
[x_lb f_lb exitflag] = quadprog(2*newH, newc, Aineq, bineq, newA, b, ...
    newLB, newUB, [], opt);

%% Computing the upper bound of the original problem with constrained minimization
% Objective function for the constrained minimization
obj_func = @(x) (x'*H*x + c'*x);
opt = optimset('Display','off');
% Nonlinear constrains must be given by a user-defined function nonlin_con
[x_ub f_ub exitflag] = fmincon(obj_func, zeros(size(ub)), [], [], [], [], ...
    lb, ub, @(x)nonlin_con(x,Qm, A, b), opt);

%% Calculate the boundaries for the auxiliary variables
function [vub, vlb] = calc_bounds(ub,lb,i,j)
vub = ub(i)*ub(j);
vlb = lb(i)*lb(j);

%% Calculate the nonlinear constrains
function [c, ceq] = nonlin_con(x,Qm,A,b)
% No nonlinear inequalities, therefore output c = 0
c = 0;
% Calculate the nonlinear equality constrains
ceq = zeros(size(Qm));
% f_ineq(x) = x'*Qi*x + ai*x - bi
for i = 1:numel(ceq)
    ceq(i) = x'*Qm{i}*x + A(i,:)*x - b(i);
end


