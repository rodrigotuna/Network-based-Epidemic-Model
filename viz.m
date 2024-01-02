%% Plots and Visualizations of Results 
% Author: Rodrigo Tuna
% e-mal: rodrigo.tunade@studio.unibo.it
% Matricola: 190011317

T = 300;
% Topology analysis
% Erdos-Renyi
figure(1);
for p = [0.10, 0.25, 0.50, 0.75]
    mat = sim_mat("erdos", 1.3, 0.999, T, 0.01, "random", "none", 12008, p);
    plot_state(mat.sim_matrix, 1, 150);
    hold on
end
hold off
legend("p=0.10","p=0.25", "p=0.50", "p=0.75");
exportgraphics(figure(1), 'plots/erdos.pdf');
%Watts-Strogatz
figure(2);
for beta = [0.20, 0.40, 0.60]
    mat = sim_mat("watts", 1.3, 0.999, T, 0.01, "random", "none", 12008, 25, beta);
    plot_state(mat.sim_matrix, 1, 150);
    hold on
end
hold off
legend("beta = 0.20", "beta = 0.40", "beta = 0.60");
exportgraphics(figure(2), 'plots/watts1.pdf');


figure(2);
for k = [2, 5, 10, 25, 50]
    mat = sim_mat("watts", 1.3, 0.999, T, 0.01, "random", "none", 12008, k, 0.20);
    plot_state(mat.sim_matrix, 1, 300);
    hold on
end
hold off
legend("k = 2", "k = 5", "k = 10", "k = 25","k = 50");
exportgraphics(figure(2), 'plots/watts2.pdf');

%All
figure(3);
mat = sim_mat("cite", 1.3, 0.999, T, 0.01, "random", "none");
plot_state(mat.sim_matrix, 1, 200);
hold on
mat = sim_mat("asia", 1.3, 0.999, T, 0.01, "random", "none");
plot_state(mat.sim_matrix, 1, 200);
hold on 
mat = sim_mat("erdos", 1.3, 0.999, T, 0.01, "random", "none", 12008, 0.50);
plot_state(mat.sim_matrix, 1, 200);
hold on 
mat = sim_mat("watts", 1.3, 0.999, T, 0.01, "random", "none", 12008, 25, 0.20);
plot_state(mat.sim_matrix, 1, 200);
hold on  
mat = sim_mat("watts", 1.3, 0.999, T, 0.01, "random", "none", 12008, 5, 0.20);
plot_state(mat.sim_matrix, 1, 200);
hold off
legend("collaboration", "LastFM", "erdos p=0.50","watts k=25", "watts k=5");
exportgraphics(figure(3), 'plots/all.pdf');

% Patient zero analysis
for data = ["cite", "asia"]
    %Impact of the number of infected
    figure(4);
    for p = [0.001, 0.01, 0.1]
        mat = sim_mat(data, 1.3, 0.999, T, p, "random", "none");
        plot_state(mat.sim_matrix, 1, 150);
        hold on
    end
    hold off
    legend("p0 = 0.001", "p0 = 0.01", "p0 = 0.1");
    exportgraphics(figure(4), 'plots/p0_' + data +'.pdf');

    figure(5);
    for strat = ["random", "best"]
        mat = sim_mat(data, 1.3, 0.999, T, 0.001, strat, "none");
        plot_state(mat.sim_matrix, 1, 75);
        hold on
    end

    hold off
    legend("random", "best");
    exportgraphics(figure(5), 'plots/strat_' + data +'.pdf');
    figure(6);
    % Parameter variation
    for im_decay = [1,0.9999, 0.999, 0.99, 0.9]
        mat = sim_mat(data, 1.3, im_decay, T, 0.01, "random", "none");
        plot_state(mat.sim_matrix, 1, 150);
        hold on
    end
    hold off
    legend("decay=1", "decay=0.9999", "decay=0.999", "decay=0.99", "decay=0.9");
    exportgraphics(figure(6), 'plots/im_decay_' + data +'.pdf')
    
    figure(7);
    for recov_rate = [1.1, 1.3, 1.5, 1.8]
        mat = sim_mat(data, recov_rate, 0.999, T, 0.01, "random", "none");
        plot_state(mat.sim_matrix, 1, 150);
        hold on
    end
    hold off
    legend("rate = 1.1", "rate = 1.3", "rate = 1.5", "rate = 1.8");
    exportgraphics(figure(7), 'plots/rec_rate_' + data +'.pdf');
    
    figure(8);
    % Efficacy of containment measures 
    for cont = ["none", "mask", "quarantine", "vaccine"]
        mat = sim_mat(data, 1.3, 0.999, T, 0.01, "random", cont);
        plot_state(mat.sim_matrix, 1,200);
        hold on
    end
    hold off
    legend("none", "mask", "quarantine", "vaccine");
    exportgraphics(figure(8), 'plots/measures_' + data +'.pdf');

    figure(9);
    for st = [0,1,2]
        mat = sim_mat(data, 1.3, 0.999, T, 0.01, "random", "none");
        plot_state(mat.sim_matrix, st,301);
        hold on
    end
    hold off
    legend("susceptible", "infected", "immune");
    ylabel("Percentage of nodes");
    exportgraphics(figure(9), 'plots/all_hist2_' + data +'.pdf');
end

function mat = sim_mat(top, recover_rate, immunity_decay, T, infected_0, initial_strat,cont,varargin)
    if top == "erdos"
        top_string = top + "_" + varargin{1} + "_" + varargin{2};
    elseif top == "watts"
        top_string = top + "_" + varargin{1} + "_" + varargin{2} + varargin{3};
    else
        top_string = top;
    end
    file_str = "output/" + top_string + '_' + recover_rate + '_' +  immunity_decay + '_' + ...
        initial_strat + infected_0 + '_' + cont + '_simulation.mat';
    mat = load(file_str);
end

function plot_state(matrix, state, until) 
    state = matrix == state;
    y = sum(state, 1) / size(matrix, 1);
    plot(y(1:until));
    xlabel("Iteration");
    ylabel("Percentage of infected nodes");
end