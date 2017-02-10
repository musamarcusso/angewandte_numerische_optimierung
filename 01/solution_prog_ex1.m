clc
close all
clear all

data2;
[mass, costs] = coal_problem(Pel, ns, np, etap, Hs, qs, asmax);

figure('windowstyle', 'docked')
bar(mass')
labels = cell(1,ns);
for i = 1:ns
    	labels{i} = sprintf('Coal mine #%d', i);
end
legend(labels)
grid on

title(sprintf('Results of the linear optimization problem (Total cost: $%4.2f)', costs), ...
    'FontSize', 12)
xlabel('Power plant index', 'FontSize', 12);
ylabel('Amount of coal [ton]', 'FontSize', 12)