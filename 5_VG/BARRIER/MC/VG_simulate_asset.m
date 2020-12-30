function [S,ST]=VG_simulate_asset(par,Nsim,M)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   INPUT: - par = parameters structure{
    %                  - par.sigma 
    %                  - par.theta 
    %                  - par.kVG
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
    kVG = par.kVG;
    dt = T/M;
    a = dt/kVG;
    
    %% Compute drift in Q-dynamics
    psi = @(u) -1/kVG .* log(1+ 0.5*(u.^2.*sigma.^2.*kVG)-1i*theta.*kVG.*u);
    drift=r-psi(-1i); % risk neutral drift
   
    %% Simulation
    X=zeros(Nsim,M+1); Z=randn(Nsim,M);
    for i=1:M 
        % Sample the Gamma subordinator
        dSub = kVG*icdf('Gamma',rand(Nsim,1),a,1);
        % Sample the process
        X(:,i+1) = X(:,i)+drift*dt+theta*dSub+sigma*sqrt(dSub).*randn(Nsim,1);
    end
       
    S  = S0*exp(X);
    ST = S(:,end);
end          