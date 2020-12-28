%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put Options in under Kou model using 
% Carr and Madan Formula via FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
% Market parameters
r = 0.01;               % riskfree interest rate 
S0 = 100;            % spot price
% Model parameters
sigma = 0.3;         % standard deviation 
p = 0.5;             % pos jump prob
lambdaK = 2;         % jump time intensity
lambdap = 15;
lambdam = 18;
% Contract parameters
T = 1;                  % maturity
Strike = [S0; S0*0.9];                 % strike

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'p',p,'lambdap',lambdap,'lambdam',lambdam,'lambdaK',lambdaK);

% Discretization parameter
Npow = 20; A = 1200;

% Chararacteristic function
CharFunc = @(v) exp(T* char_exponent_Kou(v,par) );

% Carr&Madan method
disp("Vanilla prices using Carr-Madan under Kou model:")
[call_price, put_price] = fun_FFT_CM(par, Strike, CharFunc, Npow,A)
