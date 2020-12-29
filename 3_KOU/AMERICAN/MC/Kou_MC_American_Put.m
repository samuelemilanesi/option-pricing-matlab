%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price an American put option under Kou model 
% using plain MC and Longstaff-Schwartz projection method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Parameters
% Market parameters
r = 0.001;              % riskfree interest rate 
S0 = 1;                 % spot price
% Model parameters
sigma = 0.4;            % BM vol
p = 0.5;                % jump vol
lambdap = 5;            % pos jump intensity
lambdam = 6;            % neg jump intensity
lambdaK = 3;            % n jumps intensity
% Contract parameters
T = 1;                  % maturity
K = S0;                 % strike
M = round(12*T);        % monthly monitoring
dt= T/M;
disc_payoff_fun = @(S) exp(-r*T).*max(mean(S,2)-K,0);     % disc payoff function. S(i,:) = i-th simulation of an underlying PATH 

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'p',p,'lambdap',lambdap,'lambdam',lambdam,'lambdaK',lambdaK);

% Discretization parameter
Nsim = 1e6;             % number of MC simulations 
%% Underlying asset simulation
S = Kou_simulate_asset(par, Nsim,M);
%% Longstaff-Schwartz
ExerciseTime = M*ones(Nsim,1);
Cashflow = disc_payoff_fun(S);
for step=M-1:-1:1
    InMoney=find( S(:,step+1)<K );              % select only in-the-money simultion at time step+1
    Stemp=S(InMoney,step+1);
    %-------- Regression on the span{1, S, S.^2}
    RegressionMatrix=[ones(size(Stemp)), Stemp, Stemp.^2];
    YData=Cashflow(InMoney).*exp( -r*dt*(ExerciseTime(InMoney)-step));
    alpha=RegressionMatrix\YData;
    %----------------------------------------------------------------------
    CV=RegressionMatrix*alpha;                  % Continuation value: value if I hold the option
    IV=K-Stemp;                                 % Intrinsic value: value if I exercise the option
    %-------- In the case of Early Exercise
    EarlyExercise_inMoney_index=find( IV> CV );                     % select the simulations in which I exercise the option at time `step`
    EarlyExercise_index=InMoney(EarlyExercise_inMoney_index);       % among those, consider only the one in the money
    Cashflow(EarlyExercise_index)=IV(EarlyExercise_inMoney_index);  % compute the cashflow of the exercise (i.e. the IV)
    ExerciseTime(EarlyExercise_index)=step;                         % save that the exercise time is `step`
end
[Kou_american_put_price,~,CI]=normfit( Cashflow.*exp(-r*dt*ExerciseTime) )