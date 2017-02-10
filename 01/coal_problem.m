%% coal_problem 
% Function to generate the matrices to calculate the distribution of coal
% between np power plants from ns coal mines using linear programming
% 
% Function arguments
% Pel       Minimum desired overall eletric power to be produced [Joule]
% ns         Number of coal mines available 
% np        Number of power plants available
% etap    Power plants efficiency vector [%]
% Hs        Vector of heating values for each mine [Joule/ton]
% qs         Vector of coal prices [$/ton]
% asmax Vector of maximum coal production for each mine [ton]
function [mass, costs] = coal_problem(Pel, ns, np, etap, Hs, qs, asmax)

% Verify if the minimum desired energy can be generated from all
% the coal available

max_efficiency = max(etap(:));
max_energy = max_efficiency * asmax' * Hs;

mass = [];
costs = [];

if max_energy < Pel
        error('Error - The required amount of energy cannot be produced')
        return;
end
disp(sprintf('Number of mines: %d', ns))
disp(sprintf('Number of power plants: %d', np))
disp(sprintf('Degrees of freedom: %d', ns*np - (ns +1)))
disp(' ')

% Generate the coefficients of the linear objective function c
c = [];
for i = 1:ns
    c = vertcat(c, qs(i) * ones(np,1));
end
c = double(c);
disp('Transposed vector of coefficients of the linear objective function')
disp('Coal price qs')
display(c')

% Generate the inequality contrains matrices (A*x <= b)
A = generateA(ns,np, Hs, etap);
b = generateB(ns, Pel, asmax);

format short
format compact
disp(' ')
disp('Inequality constrains')
disp('Line one: minimum amount of energy that should be produced')
disp('Other lines: maximum amount of coal that can be delivered by each coal mine')
display(A)
display(b)

% Convert to double precision 
A = double(A);
b = double(b);

% Set the option to display the results at each iteration
optset = optimset('Display', 'iter');
% Run linear programming routine with the addtional constraint that the
% masses should always be positive
x_final = linprog(c, A, b, [], [], zeros(size(c)), [], [], optset);

% Generate the output matrices 
%  Mass matrix: the calculated amount of coal that should be delivered at
%  each power plant
mass = zeros(ns,np);
count = 1;
for i = 1:ns
    for j = 1:np
        mass(i,j) = x_final(count);
        count = count + 1;
    end
end

% Costs: total cost for producing the energy
costs = 0;
for i = 1:ns
    costs = costs + sum(mass(i,:) * qs(i));
end

% Calculation the energy to be produced
energy = mass;
for i = 1:ns
        energy(i,:) = energy(i,:) * Hs(i);
end
for i = 1:np
        energy(:,i) = energy(:,i) * etap(i);
end

disp('Masses consumed by the power plants')
disp(mass)

disp('Total cost')
disp(costs)

disp('Total energy provided')
disp(sum(energy(:)))

%% generateA
% Generates the inequality matrix A for the linear programming routine.
% The signal of the inequality constraint is the opposite of the adopted
% standard presented on lectures, here the linear programming routine
% requires A * x <= b
function outA = generateA(ns,np, Hs, etap)
% Calculation of the first line (el) of the matrix to satisfy the minimum amount
% of energy that should be produced
outA = zeros(1+ns, ns*np);
el = Hs*etap';
outA(1,:) = -1 * reshape(el', 1, ns*np);
% Generation of the rest of the matrix's lines. Each line will represent
% the constraint of maximum amount of coal that can be produced at each
% coal mine
aux = zeros(ns, np);
aux_mat = [];
aux(1,:) = 1;
for i = 1:ns
    if isempty(aux_mat)
        aux_mat = circshift(aux, i - 1);
    else
        aux_mat = horzcat(aux_mat, circshift(aux, i - 1));
    end
end

outA(2:end, :) = aux_mat;
return

%% generateB
% Generation of the column vector b of the inequality constraints
function outB = generateB(ns, Pel, asmax)
outB = zeros(1+ns, 1);
outB(1) = -Pel;
outB(2:end,1) = asmax';
return