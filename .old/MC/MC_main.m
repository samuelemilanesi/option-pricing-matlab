clear all
close all
clc

%% INPUT PARAMETERS 
% Market parameteres
r = 0.2
S0 = 100

% Contract parameters
T = 1 
Stike = 100

% Model parameters
sigma = 0.4

% Numerical parameters 
Nsim = 1e7 

%% BUILDING BLOCKS SIMULATION 
Z = randn(Nsim,1); 


%% MODEL FUNCTIONS 
simulate_ST = @(n, seed) simulate_ST_BS(r, S0, T, Strike, sigma, n, seed);
simulate_prices = @(ST) prices_EU_CALL(ST,Strike);

% simulate_prices_AV = @d


%% UNDERLYING SIMULATION  
ST = simulate_ST(Nsim,0);

%% SIMULATE DISCOUNTED PAYOFF
simPrice = simulate_prices(ST)

%% Price evaluation 
Price = normfit(discounted_payoff) 

