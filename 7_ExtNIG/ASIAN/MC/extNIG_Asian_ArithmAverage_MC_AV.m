%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price Asian Arithmetic Average Call option under NIG model 
% using MC with Antithetic Variance Reduction
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

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'sigmaGBM',sigmaGBM,'theta',theta,'kNIG',kNIG);
% Discretization parameter
Nsim = 1e5;             % number of MC simulations 

%% Simulate Underlying Asset
% Simulate the underlying at time T
[S, SAV] = extNIG_simulate_assetAV(par,Nsim,M);
%% Compute the discounted payoff
DiscPayoff = disc_payoff_fun(S);
DiscPayoffAV = disc_payoff_fun(SAV);

%% Compute call price and asymptotic CI 
disp("extNIG Model - Asian Arithmetic Average call option price MC with Antithetic Variance Reduction:")
[asian_call_price_AV, ~, asian_call_CI_price_AV] = normfit((DiscPayoff+DiscPayoffAV)/2)

CI_len = asian_call_CI_price_AV(2)-asian_call_CI_price_AV(1)