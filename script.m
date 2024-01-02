%% Runs all experiments
% Author: Rodrigo Tuna
% e-mal: rodrigo.tunade@studio.unibo.it
% Matricola: 190011317

T = 300;
% Topology analysis

% Erdos-Renyi
for p = [0.10, 0.25, 0.50, 0.75]
    simulation("erdos", 1.3, 0.999, T, 0.01, "random", "none", 12008, p);
end

%Watts-Strogatz
for beta = [0.20, 0.40, 0.60]
    simulation("watts", 1.3, 0.999, T, 0.01, "random", "none", 12008, 25, beta);
end

for k = [2, 5, 10, 50]
    simulation("watts", 1.3, 0.999, T, 0.01, "random", "none", 12008, k, 0.20);
end

%Real World Data 
simulation("cite", 1.3, 0.999, T, 0.01, "random", "none");
simulation("asia", 1.3, 0.999, T, 0.01, "random", "none");

disp("Top over")
% Patient zero analysis
for p = [0.001, 0.1]
    simulation("cite", 1.3, 0.999, T, p, "random", "none");
    simulation("asia", 1.3, 0.999, T, p, "random", "none");
end

for p = [0.001, 0.01, 0.1]
    simulation("cite", 1.3, 0.999, T, p, "best", "none");
    simulation("asia", 1.3, 0.999, T, p, "best", "none");
end

disp("patient over")
% Parameter variation
for im_decay = [1,0.9999, 0.99, 0.9]
    simulation("cite", 1.3, im_decay, T, 0.01, "random", "none");
    simulation("asia", 1.3, im_decay, T, 0.01, "random", "none");
end

for recov_rate = [1.1, 1.5, 1.8]
    simulation("cite", recov_rate, 0.999, T, 0.01, "random", "none");
    simulation("asia", recov_rate, 0.999, T, 0.01, "random", "none");
end

disp("param over")

% Efficacy of containment measures
for cont = ["mask", "quarantine", "vaccine"]
    simulation("cite", 1.3, 0.999, T, 0.01, "random", cont);
    simulation("asia", 1.3, 0.999, T, 0.01, "random", cont);
end