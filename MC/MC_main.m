%% INPUT PARAMETERS 
% Market parameteres
r = 0.2
S0 = 100

% Contract parameters
T = 1 
Stike = 100

% Model parameters
sigma = 0.4

% Discretization parameters 
Nsim = 1e7 

%%  Undelying simulation 

ST = simulate_ST_BS(r, S0, T, Strike, sigma, Nsim, seed=0);% if path independent
%S = simulate pathblabla % if path dependent 

%% Simulate discounted payoff 

discounted_payoff = exp(-r*T)*payoff_EU_Call(ST, Strike)

%% Price evaluation 

Price = normfit(discounted_payoff) 

