%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price Asian Arithmetic Average Call option under NIG model using MC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Parameters
% Market parameters
r = 0.01;            % riskfree interest rate 
S0 = 100;            % spot price
% Model parameters
sigma = 0.6;
sigmaGBM=0.4;
theta = 0.05;
kNIG   = 0.2;
% Contract parameters
T = 1;                  % maturity
K = S0;                 % strike
M = round(12*T);        % monthly monitoring
disc_payoff_fun = @(S) exp(-r*T).*max(mean(S,2)-K,0);     % disc payoff function. S(i,:) = i-th simulation of an underlying PATH 

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'sigmaGBM',sigmaGBMtheta',theta,'kNIG',kNIG);
% Discretization parameter
Nsim = 1e6;             % number of MC simulations 

%% Simulate Underlying Asset
% Simulate the underlying at time T
S = extNIG_simulate_asset(par,Nsim,M);
%% Compute the discounted payoff
DiscPayoff = disc_payoff_fun(S);

%% Compute call price and asymptotic CI 
disp("NIG Model - Asian Arithmetic Average call option price via plain MC:")
[asian_call_price, ~, asian_call_CI_price] = normfit(DiscPayoff)

CI_len = asian_call_CI_price(2)-asian_call_CI_price(1)