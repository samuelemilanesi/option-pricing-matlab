%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price Lookback fixed strike call option under VG model using MC
% with Antithetic Variance Reduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Parameters
% Market parameters
r = 0.01;            % riskfree interest rate 
S0 = 100;            % spot price
% Model parameters
sigma = 0.6;
theta = 0.05;
kVG   = 0.2;
% Contract parameters
T = 1;                  % maturity
K = S0;                 % strike
M = round(12*T);        % monthly monitoring

disc_payoff_fun = @(S) exp(-r*T).*max(max(S,[],2)-K,0);     % disc payoff function. S(i,:) = i-th simulation of an underlying PATH 
par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'theta',theta,'kVG',kVG);
% Discretization parameter
Nsim = 1e6;             % number of MC simulations 

%% Simulate Underlying Asset
% Simulate the underlying at time T
[S, SAV] = VG_simulate_asset_AV(par,Nsim,M);
%% Compute the discounted payoff
DiscPayoff = disc_payoff_fun(S);
DiscPayoffAV = disc_payoff_fun(SAV);

%% Compute call price and asymptotic CI 
disp("VG Model - Lookback fixed strike call option price MC with Antithetic Variance Reduction:")
[lb_fixStrike_call_price_AV, ~, lb_fixStrike_call_price_CI_AV] = normfit((DiscPayoff+DiscPayoffAV)/2)

CI_len = lb_fixStrike_call_price_CI_AV(2)-lb_fixStrike_call_price_CI_AV(1)