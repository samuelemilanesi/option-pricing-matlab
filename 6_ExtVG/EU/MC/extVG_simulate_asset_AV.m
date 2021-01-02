function [S,SAV,ST,STAV]=extVG_simulate_asset_AV(par,Nsim,M)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   INPUT: - par = parameters structure{
    %                  - par.sigma = conditional vol
    %                  - par.sigmaGBM = uncond GBM vol 
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
    sigmaGBM = par.sigmaGBM;
    theta = par.theta;
    kVG = par.kVG;
    dt = T/M;
    a = dt/kVG;
    
    %% Compute drift in Q-dynamics
    psi = @(u) -u.^2*sigmaGBM^2/2 -log(1+ 0.5*(u.^2.*sigma.^2.*kVG)-1i*theta.*kVG.*u)./kVG;
    drift=r-psi(-1i); % risk neutral drift
   
    %% Simulation
    X=zeros(Nsim,M+1); Z=randn(Nsim,M);
    XAV = X;
    for i=1:M 
        % Sample the Gamma subordinator
        dSub = kVG*icdf('Gamma',rand(Nsim,1),a,1);
        % Sample Gaussian extension
        Zext = randn(Nsim,1);
        % Sample Gaussian subordinated
        Z = randn(Nsim,1);
        % Sample the process
        X(:,i+1) = X(:,i)+drift*dt+theta*dSub+sigma*sqrt(dSub).*Z + sigmaGBM*sqrt(dt)*Zext;
        XAV(:,i+1) = XAV(:,i)+drift*dt+theta*dSub-sigma*sqrt(dSub).*Z - sigmaGBM*sqrt(dt)*Zext;
        
    end
       
    S  = S0*exp(X);
    SAV = S0*exp(XAV);

    STAV = SAV(:,end);
    ST = S(:,end);
end          