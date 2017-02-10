clc
clear all
close all
warning off

data_prog2;

[f_lb, f_ub] = convex_bound(n, m, c, H, Q, A, b, lb, ub)