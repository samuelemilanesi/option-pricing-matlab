%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price an Asian Arithmetic Average Call option under B&S model 
% using plain MC 
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
M = round(12*T);          % monthly monitoring
disc_payoff_fun = @(S) exp(-r*T).*max(mean(S,2)-K,0);     % disc payoff function. S(i,:) = i-th simulation of an underlying PATH 

% Discretization parameter
Nsim = 1e7;             % number of MC simulations 

%% Simulate Underlying Asset
% Simulate the underlying at time T
S = BS_simulate(S0, r, sigma, T, Nsim, M);   % simulates Nsim path of M steps in time of the underlying asset

%% Compute the discounted payoff
DiscPayoff = disc_payoff_fun(S);

%% Compute call price and asymptotic CI 
disp("Asian Ar. Ave. call option price via plain MC:")
[asian_call_price, ~,asian_call_price_CI] = normfit(DiscPayoff)