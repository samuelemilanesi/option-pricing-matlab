%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price a Lookback fixed strike call option under Merton model 
% using plain MC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Parameters
% Market parameters
r = 0.0367;             % riskfree interest rate 
S0 = 100;               % spot price
% Model parameters
sigma = 0.17801;        % BM vol
delta = sqrt(0.4);      % jump vol
mu    = 0.01;           % jump drift
lambdaK = 0.2;          % n jumps intensity
% Contract parameters
T = 1;                  % maturity
K = S0;                 % strike
D = 90;                 % barrier
M = 52*T;               % weekly monitoring
disc_payoff_fun = @(S) exp(-r*T).*max(max(S,[],2)-K,0);     % disc payoff function. S(i,:) = i-th simulation of an underlying PATH 

% Discretization parameter
Nsim = 1e5;             % number of MC simulations 

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'mu',mu,'delta',delta,'lambdaK',lambdaK);

%% Simulate Underlying Asset
% Simulate the underlying at time T
S = Merton_simulate_asset(par,Nsim,M);
%% Compute the discounted payoff
DiscPayoff = disc_payoff_fun(S);

%% Compute call price and asymptotic CI 
disp("Merton model - Lookback fixed strike call option price via plain MC:")
[lb_fixedStrike_price, ~, lb_fixedStrike_price_CI] = normfit(DiscPayoff)