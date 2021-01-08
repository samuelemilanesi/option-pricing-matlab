%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price Barrier D&O Call option under NIG model using MC
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
D = S0*0.9;             % barrier
M = round(12*T);        % monthly monitoring

disc_payoff_fun = @(S) exp(-r * T) * max(S(:,end) - K, 0).*(min(S,[],2)>D);
par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'sigmaGBM,sigmaGBM,'theta',theta,'kNIG',kNIG);
% Discretization parameter
Nsim = 1e6;             % number of MC simulations 

%% Simulate Underlying Asset
% Simulate the underlying at time T
S = extNIG_simulate_assetAV(par,Nsim,M);
%% Compute the discounted payoff
DiscPayoff = disc_payoff_fun(S);

%% Compute call price and asymptotic CI 
disp("NIG Model - Barrier D&O call option price via plain MC:")
[DO_call_price, ~, DO_call_CI_price] = normfit(DiscPayoff)

CI_len = DO_call_CI_price(2)-DO_call_CI_price(1)