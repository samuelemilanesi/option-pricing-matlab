%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put Options in under Kou model using 
% Carr and Madan Formula via FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
% Market parameters
r = 0.0367;             % riskfree interest rate 
S0 = 100;               % spot price
% Model parameters
sigma = 0.17801;        % BM vol
p = 0.4;                % prob of positive jump
lambdap = 4;            % pos jump intensity
lambdam = 2;            % neg jump intensity
lambdaK = 0.2;          % n jumps intensity
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
