function [S,ST,V]=Bates_simulate_asset(par,Nsim,M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: 
% par{ epsilon  = vol-of-vol
%       k       = mean reversion speed
%       rho     = correlation
%       theta   = long-term mean
%       V0      = initial var-of-vol
%       kZ      = log(1+kZ) = mean of jump size
%       lambda  = intensity of jump time
%       sigma   = vol of jump size    }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% de-struct params 
    V0 = par.V0;
    epsilon=par.epsilon;
    S0 = par.S0;
    r = par.r;
    T = par.TTM;
    sigma = par.sigma;
    theta=par.theta;
    psi=par.kappa;
    kappa=par.kappaZ;		 
    lambda  = par.lambda;
    rho=par.rho;
    %% Compute drift in Q-dynamics
    dt=T/M;
   
    %% Simulation
    X = zeros(Nsim,M+1); X(:,1) = log(S0);
    V = zeros(Nsim,M+1); V(:,1) = V0;

    for t=1:M
        Njumps = icdf('Poisson',rand(Nsim,1),lambda*dt);
        R = [1 rho; rho 1];
        XY= mvnrnd([0;0],R,Nsim);
        x = XY(:,1);
        y = XY(:,2);
        V(:,t+1)=V(:,t)+psi*(theta-max(V(:,t),0))*dt+epsilon*sqrt(dt*max(V(:,t),0)).*y;
        X(:,t+1)=X(:,t)+(r-lambda*kappa-1/2*V(:,t))*dt+sqrt(dt*max(V(:,t),0)).*x;
        for i=1:Nsim
            for j=1:Njumps(i)
                X(i,t+1) = X(i,t+1)+ randn*sigma+log(1+kappa);
            end
        end
    end
    S=exp(X);
    ST = S(:,end);
end          