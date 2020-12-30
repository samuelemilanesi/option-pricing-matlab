%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put Options under VG model using 
% MC with Control Variable
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
K = S0;  % strike

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'theta',theta,'kVG',kVG);
% Discretization parameter
Nsim = 1e7;             % number of MC simulations 

%% Define control variable 
ExpValue_g_CV = S0*exp(r*T);      % use as control variable ST. g(ST) = ST. 
% Estimate alpha
Nsim2=1e6;                      % number of simulations used for alpha estimation
[~,g_CV] = VG_simulate_asset(par,Nsim2,1);    % g(ST) = ST
f=exp(-r*T)*max(g_CV-K,0);      % f = dis = Merton_simulate_asset(par,Nsim,1)counted payoff function 
VC=cov(f,g_CV);
alpha=-VC(1,2)/VC(2,2)          % optimal alpha 

%% Simulate Underlying Asset
[~,ST] = VG_simulate_asset(par,Nsim,1);  % simulate ST with a simul indep from the one used to estimate alpha

%% Compute the discounted payoff independetly from the simulation of f used to estimate alpha 
DiscPayoff = exp(-r * T) * max(ST - K, 0);

%% Compute call price and asymptotic CI using the control variable
disp("VG model - Vanilla options price via MC with Control Variable:")
[call_price, ~, call_CI_price] = normfit(DiscPayoff + alpha*(ST - ExpValue_g_CV))
CI_len = call_CI_price(2)-call_CI_price(1)
%% Compute put price and asymptotic CI via Call-Put Parity
put_parity = @(call_p) call_p - S0 + K*exp(-r*T);

put_price = put_parity(call_price)
put_CI_price = put_parity(call_CI_price)
