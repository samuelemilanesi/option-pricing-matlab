function [S,SAV,ST,STAV]=Merton_simulate_asset_AV(par,Nsim,M)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   MC Merton simulation using antithetic variance reduction
    %   INPUT: - par = parameters structure{
    %                  - par.S0 = underlying  spot price
    %                  - par.r = risk free rate
    %                  - par.TTM = time to maturity
    %                  - par.sigma = BM vol 
    %                  - par.mu = jump drift
    %                  - par.delta = jump vol
    %                  - par.lambdaK = jump time intensity
    %                  }
    %           - Nsim = number of simulations
    %           - M = number of time steps
    %   OUTPUT: - S = matrix of Nsim underlying paths (one for each row)
    %           - ST = vector of Nsim underlying sim at time T 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% de-struct params 
    S0 = par.S0;
    r = par.r;
    T = par.TTM;
    sigma = par.sigma;
    mu = par.mu;
    delta = par.delta;
    lambdaK  = par.lambdaK;
    
    %% Compute drift in Q-dynamics
    dt=T/M;
    psi = @(v)-sigma^2/2*v.^2+lambdaK.*(exp(-delta^2/2*v.^2+mu*1i*v)-1); % without drift   
    drift=r-psi(-1i); % risk neutral drift
   
    %% Simulation
    X = zeros(Nsim,M+1); 
    XAV = X; 
    
    for t=1:M
        Njumps = icdf('Poisson',rand(Nsim,1),lambdaK*dt);
        for i=1:Nsim
            % sample continuous part 
            Z = randn;
            X(i,t+1) = X(i,t)+drift*dt+sigma*sqrt(dt)*Z;
            XAV(i,t+1) = XAV(i,t)+drift*dt-sigma*sqrt(dt)*Z;
            % add the jump component
            for j=1:Njumps(i)
                Zjump = randn;
                X(i,t+1) = X(i,t+1)+ delta*Zjump + mu;
                XAV(i,t+1) = XAV(i,t+1)+ delta*Zjump + mu;
            end
        end
    end
    S = S0*exp(X);
    ST = S(:,end);
    SAV = S0*exp(XAV);
    STAV = SAV(:,end);
end          