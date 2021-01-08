%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price an American put option under extNIG model 
% using plain MC and Longstaff-Schwartz projection method
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
M = round(12*T);        % monthly monitoring
dt= T/M;
disc_payoff_fun= @(S) max(K-S(:,end),0);     % disc payoff function. S(i,:) = i-th simulation of an underlying PATH 
%disc_payoff_fun= @(S) max(S(:,end)-K,0);
par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'sigmaGBM',sigmaGBM,'theta',theta,'kNIG',kNIG);

% Discretization parameter
Nsim = 1e5;             % number of MC simulations 
%% Underlying asset simulation
S = extNIG_simulate_asset(par, Nsim,M);
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
[extNIG_american_put_price,~,CI]=normfit( Cashflow.*exp(-r*dt*ExerciseTime) )

CI_len = CI(2)-CI(1)