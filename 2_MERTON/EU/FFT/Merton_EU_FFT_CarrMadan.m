%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put Options in under Merton model using 
% Carr and Madan Formula via FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
Strike = S0;                 % strike

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'mu',mu,'delta',delta,'lambdaK',lambdaK);

% Discretization parameter
Npow = 20; A = 1200;

% Chararacteristic function
CharFunc = @(v) exp(T* char_exponent_Merton(v,par) );

% Carr&Madan method
disp("Vanilla prices using Carr-Madan under Merton model:")
[call_price, put_price] = fun_FFT_CM(par, Strike, CharFunc, Npow,A)
