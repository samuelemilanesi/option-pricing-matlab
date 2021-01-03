clear; close all;

%% Parameters
S_0 = 100;
% Market parameters
param.rf = 0.0367;             % riskfree interest rate 
% Model parameters
param.sigma = 0.17801;        % BM vol
param.delta = sqrt(0.4);      % jump vol
param.mu    = 0.01;           % jump drift
param.lambdaK = 0.2;          % n jumps intensity
% Contract parameters
param.T = 1;                  % maturity
K = S_0;                 % strike
M = round(52*param.T);        % weekly monitoring
param.q= 0; %dividend
param.dt= param.T/M;

param.distr=2; % merton
Barrier=0.9*S_0; N=2^12;
[S,v] = CONV( S_0, K, M, N, Barrier, param);
price=interp1(S,v,S_0,'spline')

MC_gap = price - 10.369088217324160 
close all