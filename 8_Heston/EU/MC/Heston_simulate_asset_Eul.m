function [S,ST]=Heston_simulate_asset_Eul(par,Nsim,M)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   INPUT: - par = parameters structure{
    %                  - par.S0 = underlying  spot price
    %                  - par.r = risk free rate
    %                  - par.TTM = time to maturity
    %                  - par.epsilon = vol-of-vol
    %                  - par.kappa   = mean reversion speed
    %                  - par.rho     = correlation
    %                  - par.V0      = initial var
    %                  - par.theta   = mean
    %                  }
    %           - Nsim = number of simulations
    %           - M = number of time steps
    %   OUTPUT: - S = matrix of Nsim underlying paths (one for each row)
    %           - ST = vector of Nsim underlying sim at time T 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% de-struct params 
    epsilon = par.epsilon; % vol-of-vol
    k = par.kappa; % mean reversion speed
    rho = par.rho; % correlation
    theta = par.theta; % mean 
    V0 = par.V0;
    T = par.TTM; 
    S0 = par.S0;
    r = par.r;
    dt=T/M; 
    V = zeros(M+1,Nsim); S = zeros(M+1,Nsim);
    
    %--- DISCRETIZE V ---------------------------------------------------------
    V(1,:)=V0;
    
 S(1,:)=S0;
 V(1,:)=V0;
    for i=1:M
        R = [1 rho; rho 1];
        XY= mvnrnd([0;0],R,Nsim);
        x = XY(:,1);
        y = XY(:,2)
        V(i+1,:)=V(i,:)+k*(theta-max(V(i,:),0))*dt+epsilon*sqrt(dt*max(V(i,:),0)).*y';
        S(i+1,:)=exp(log(S(i,:))-0.5*max(V(i,:),0)*dt+sqrt(dt*max(V(i,:),0)).*x');
    end

    S = S';
    ST = S(:,end);
end

    