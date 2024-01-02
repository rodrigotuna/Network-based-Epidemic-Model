%% Implementation of the Simulation
% Author: Rodrigo Tuna
% e-mal: rodrigo.tunade@studio.unibo.it
% Matricola: 190011317

function [] = simulation(top, recover_rate, immunity_decay, T, infected_0, initial_strat,cont,varargin)
    epsilon = 0.01; %recover probability after infection
    threshold = 0.3;
    
    
    %Topology 
    if top == "erdos"
        Adj = erdosRenyi(varargin{1},varargin{2});
        top_string = top + "_" + varargin{1} + "_" + varargin{2};
    elseif top == "watts"
        Adj = wattsStrogatz(varargin{1},varargin{2}, varargin{3});
        top_string = top + "_" + varargin{1} + "_" + varargin{2} + varargin{3};
    elseif top == "cite"
        data = readmatrix("data/ca-HepPH_index1.csv");
        G = graph(data(:,1),data(:,2));
        Adj = adjacency(G);
        top_string = top;
    elseif top == "asia"
        data = readmatrix("data/lastfm_asia_edges.csv");
        G = graph(data(:,1) + 1,data(:,2) + 1);
        Adj = adjacency(G);
        top_string = top;
    end

    N = length(Adj);
 
    %Node state 0-healthy, 1-infected, 2-immune
    state = zeros(N,1);
    
    %Node recover probability;
    recover = zeros(N,1);

    %Immunity of each node, inittialy no node is imune to the disease
    immunity = zeros(N,1);

    %Node transmissibily 
    % Added random wheights to nodes transmissibily inicating frequency of
    % contacts
    R = rand(N);
    A = Adj .* ((R + R')/2);
    
    %Select first nodes to be infected according to 2 different strategies
    %Randomly selected
    %"Best" - most central nodes
    M = round(N * infected_0);
    if initial_strat == "best"
        G = graph(Adj);
        pgr = centrality(G,'pagerank');
        [~, isort] = sort(pgr, 'descend');
        state(isort(1:M)) = 1;
        recover(isort(1:M)) = epsilon;
    elseif initial_strat == "random"
        perm = randperm(N);
        state(perm(1:M)) = 1;
        recover(perm(1:M)) = epsilon;
    end     
    
    %matrix with the history of states for analysis 
    sim_matrix = zeros(N, T+1);
    sim_matrix(1:N, 1) = state;

    %Main simulation loop
    for t = 1:T
        recovered = rand(N,1) < recover & state==1; %Nodes recovered at this iteration 
        state(recovered) = 2;
        immunity(recovered) = 1; %full immunity after recovery
        
        recover = recover * recover_rate; %update recover period;
        recover = min(recover, 1);
        
        %nodes become susceptible again
        susceptible = rand(N,1) > immunity & state == 2;
        immunity = immunity * immunity_decay;
        state(susceptible) = 0; 
        
        %Change weights of to represent syntoms of desease
        iteration_weights = ones(N,1);
        
        recover_threshold = recover-threshold;
        cond_ge = state == 1 & recover >= threshold;
        cond_le = state == 1 & recover < threshold;
        iteration_weights(cond_ge) = recover_threshold(cond_ge);
        iteration_weights(cond_le) = -recover_threshold(cond_le);
        A_it = iteration_weights .* A;

    
        Trans = (rand(N) < A_it) + eye(N); %randomly project matrix and add diag to maintain state
        new_state = Trans' * (state==1);
        new_state = min(new_state, 1);
        
        %set the recover probability of 
        recover(new_state > state) = epsilon; %small probability of recovery in next step
        new_state(state==2) = 2;
        state = new_state;

        sim_matrix(1:N, t+1) = state;
        
        %Containment measures
        if t == 5
            if cont == "mask"
                A = A * 0.5;
            elseif cont == "quarantine"
                A = (A > 0.8) .* A;
            elseif cont == "vaccine"
                M = round(N*0.5);
                perm = randperm(N);
                vaccinated = perm(1:M);
                state(vaccinated) = 2;
                immunity(vaccinated) = 1;
            end
        end
    end
    %Save states at all iterations and the graph
    save("output/" + top_string + '_' + recover_rate + '_' +  immunity_decay + '_' + ...
        initial_strat + infected_0 + '_' + cont + '_simulation.mat','sim_matrix');
    save("output/" + top_string + '_' + recover_rate + '_' +  immunity_decay + '_' + ...
        initial_strat + infected_0 + '_' + cont + '_graph.mat', 'A');
end


%Implements the Erdos-Renyi random graph model
function A = erdosRenyi(N,p)
    A = rand(N);
    A = (A + A')/2;
    A = (A < p) .* (ones(N) - eye(N));
end

%Implements the Watts-Strogatz small-world model
function A = wattsStrogatz(N,K, beta)
    A = zeros(N);
    for i = 1:N
        for k = 1:K
            j = mod(i + k, N) + 1;
            if rand(1,1) < beta
                j = randi(N-1);
                j = j + (j >= i);
            end
            A(i,j) = 1;
            A(j,i) = 1;
        end
    end
end
