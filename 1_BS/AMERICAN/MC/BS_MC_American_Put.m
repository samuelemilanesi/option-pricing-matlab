%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price an American Call option under B&S model 
% using plain MC and Longstaff-Schwartz projection method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Parameters
% Market parameters
r = 0.02;               % riskfree interest rate 
S0 = 218.75;            % spot price
% Model parameters
sigma = 0.4;            % standard deviation 
% Contract parameters
T = 1;                  % maturity
K = S0;                 % strike
D = 0.98*S0;            % barrier
M = round(12*T);        % monthly monitoring
dt= T/M;
disc_payoff_fun = @(S) exp(-r*T).*max(S(:,end)-K,0).*(min(S,[],2)>D);     % disc payoff function. S(i,:) = i-th simulation of an underlying PATH 

% Discretization parameter
Nsim = 1e6;             % number of MC simulations 

%% Underlying asset simulation
S = BS_simulate(S0,r,sigma,T,Nsim,M);
%% Longstaff-Schwartz
ExerciseTime=M*ones(Nsim,1);
Cashflow=max( K-S(:,end),0 );
for step=M-1:-1:1
    InMoney=find( S(:,step+1)<K );              % select only in-the-money simultion at time step+1
    % InMoney = find( S(:,step+1)>K );          % in case of American call option 
    Stemp=S(InMoney,step+1);
    %-------- Regression on the span{1, S, S.^2}
    RegressionMatrix=[ones(size(Stemp)), Stemp, Stemp.^2];
    YData=Cashflow(InMoney).*exp( -r*dt*(ExerciseTime(InMoney)-step));
    alpha=RegressionMatrix\YData;
    %----------------------------------------------------------------------
    CV=RegressionMatrix*alpha;                  % Continuation value: value if I hold the option
    IV=K-Stemp;                                 % Intrinsic value: value if I exercise the option
    %IV= Stemp-K;         % in case of American call option
    %-------- In the case of Early Exercise
    EarlyExercise_inMoney_index=find( IV> CV );                     % select the simulations in which I exercise the option at time `step`
    EarlyExercise_index=InMoney(EarlyExercise_inMoney_index);       % among those, consider only the one in the money
    Cashflow(EarlyExercise_index)=IV(EarlyExercise_inMoney_index);  % compute the cashflow of the exercise (i.e. the IV)
    ExerciseTime(EarlyExercise_index)=step;                         % save that the exercise time is `step`
end
[Price,~,CI]=normfit( Cashflow.*exp(-r*dt*ExerciseTime) )