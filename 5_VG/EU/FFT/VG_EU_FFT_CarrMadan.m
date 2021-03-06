%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put Options in under VG model using 
% Carr and Madan Formula via FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
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
Strike = [S0; S0*0.9];  % strike

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'theta',theta,'kVG',kVG);

% Discretization parameter
Npow = 20; A = 1200;

% Chararacteristic function
CharFunc = @(v) exp(T* char_exponent_VG(v,par) );

% Carr&Madan method
disp("Vanilla prices using Carr-Madan under VG model:")
[call_price, put_price] = fun_FFT_CM(par, Strike, CharFunc, Npow,A)
