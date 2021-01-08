function [S,ST]=NIG_simulate_asset(par,Nsim,M)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   INPUT: - par = parameters structure{
    %                  - par.sigma = conditional vol 
    %                  - par.theta = conditional part of the drift
    %                  - par.kVG   = param k
    %                  - par.S0 = underlying  spot price
    %                  - par.r = risk free rate
    %                  - par.TTM = time to maturity
    %                  }
    %           - Nsim = number of simulations
    %           - M = number of time steps
    %   OUTPUT: - S = matrix of Nsim underlying paths (one for each row)
    %           - ST = vector of Nsim underlying sim at time T 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% params 
    S0 = par.S0;
    r = par.r;
    T = par.TTM;
    sigma = par.sigma;
    theta = par.theta;
    k = par.kNIG;
    dt = T/M;
     X = zeros(Nsim, M+1);
    
   %%
   V_noDrift = @(z) (1/k)*(1 - (1 +z.^2*sigma^2*k - 2*1i*theta*k.*z).^0.5); %without drift NIG
   %V_noDrift = @(z) -log(1+z.^2*sigma^2*k/2-1i*theta*k.*z)/k; % without drift VG
   drift_rn = -V_noDrift(-1i);

   lambda=(dt^2)/k;  mu = dt;  %parametri di simulazione del subordinator
   %a=dt/k;

for j=1:M
    
    deltaS = icdf('InverseGaussian', rand(Nsim,1), mu, lambda);
    %deltaS = k*icdf('Gamma', rand(Nsim,1), a, 1);
    W = sqrt(deltaS).*randn(Nsim,1);
    
    X(:,j+1) = X(:,j) + theta*deltaS + sigma*W +(r + drift_rn)*dt;
end

       
    S  = S0*exp(X);
    ST = S(:,end);
end          