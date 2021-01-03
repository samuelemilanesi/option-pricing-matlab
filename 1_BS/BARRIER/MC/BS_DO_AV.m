%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price a Barrier Down&Out call option under B&S model 
% using MC with Antithetic Variance Reduction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Parameters
% Market parameters
r = 0.02;               % riskfree interest rate 
S0 = 218.75;            % spot price
% Model parameters
sigma = 0.2516;         % standard deviation 
% Contract parameters
T = 1;                  % maturity
K = S0;                 % strike
D = 0.98*S0;            % barrier
M = round(12*T);        % monthly monitoring
disc_payoff_fun = @(S) exp(-r*T).*max(S(:,end)-K,0).*(min(S,[],2)>D);     % disc payoff function. S(i,:) = i-th simulation of an underlying PATH 

% Discretization parameter
Nsim = 3e7;             % number of MC simulations 

%% Simulate Underlying Asset
% Simulate the underlying at time T
[S,SAV] = BS_simulate_AV(S0, r, sigma, T, Nsim, M);   % simulates Nsim path of M steps in time of the underlying asset

%% Compute the discounted payoff
DiscPayoff = disc_payoff_fun(S);            % note: handle functions are performed on columns by default, so we have to transpose the matrix
DiscPayoffAV = disc_payoff_fun(SAV);            % note: handle functions are performed on columns by default, so we have to transpose the matrix

%% Compute call price and asymptotic CI 
disp("Barrier DO option price via MC with Antithetic Variance Reduction:")
[barrier_DO_price, ~, barrier_DO_CI_price] = normfit((DiscPayoff+DiscPayoffAV)/2)