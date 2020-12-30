%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price Barrier D&O Call option under VG model using MC
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
D = S0*0.9;             % barrier
M = round(12*T);        % monthly monitoring

disc_payoff_fun = @(S) exp(-r * T) * max(S(:,end) - K, 0).*(min(S,[],2)>D);
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
disp("VG Model - Barrier D&O call option price MC with Antithetic Variance Reduction:")
[DO_call_price_AV, ~, DO_call_CI_price_AV] = normfit((DiscPayoff+DiscPayoffAV)/2)

CI_len = DO_call_CI_price_AV(2)-DO_call_CI_price_AV(1)