%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put Options in under NIG model using 
% Carr and Madan Formula via FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; 
%% Parameters
% Market parameters
r = 0.01;            % riskfree interest rate 
S0 = 100;            % spot price
% Model parameters
sigma = 0.6;
theta = 0.05;
kNIG   = 0.2;
% Contract parameters
T = 1;                  % maturity
Strike = S0;  % strike

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'theta',theta,'kNIG',kNIG);

% Discretization parameter
Npow = 20; A = 1200;

% Chararacteristic function
CharFunc = @(v) exp(T* char_exponent_NIG(v,par) );

% Carr&Madan method
disp("Vanilla prices using Carr-Madan under NIG model:")
[call_price, put_price] = fun_FFT_CM(par, Strike, CharFunc, Npow,A)
