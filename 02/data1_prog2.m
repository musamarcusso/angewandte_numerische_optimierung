c = [1 1]';
b = 0;

A = [0 0];

H = [0 0; ...
    0 1];

Q1 = [0 0.5; ...
    0.5 0];

Q = Q1;

lb = 0.1.*[1 1]';
ub = 10.* [1 1]';

n = size(H,1);
m = size(A,1);