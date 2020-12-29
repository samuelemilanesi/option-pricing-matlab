function [S,ST]=Kou_simulate_asset(par,Nsim,M)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   INPUT: - par = parameters structure{
    %                  - par.S0 = underlying  spot price
    %                  - par.r = risk free rate
    %                  - par.TTM = time to maturity
    %                  - par.sigma = BM vol 
    %                  - par.p = positive jump prob
    %                  - par.lambdap = positive jump intensity
    %                  - par.lambdam = negative jump intensity
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
    p = par.p;
    lambdap  = par.lambdap;
    lambdam  = par.lambdam;
    lambdaK  = par.lambdaK;
    
    %% Compute drift in Q-dynamics
    dt=T/M;
    psi=@(u) -sigma^2/2*u.^2+1i*u*lambdaK.*(p./(lambdap-1i*u)-(1-p)./(lambdam+1i*u));
    drift=r-psi(-1i); % risk neutral drift
   
    % Simulating -> conditional simulation
    NT=icdf('Poisson',rand(Nsim,1),lambdaK*T);
    X=zeros(Nsim,M+1); Z=randn(Nsim,M);
    for i=1:Nsim
        JumpTimes=sort(rand(NT(i),1));
        for j=1:M
            X(i,j+1)=X(i,j)+drift*dt+sigma*sqrt(dt)*Z(i,j);
            for jj=1:NT(i)
                if (JumpTimes(jj)>(j-1)*dt).*(JumpTimes(jj)<=j*dt)
                    if rand<p
                        Jump=icdf('Exponential',rand,1/lambdap);
                    else
                        Jump=-icdf('Exponential',rand,1/lambdam);
                    end
                    X(i,j+1)=X(i,j+1)+Jump;
                end
            end
        end
    end
    S=S0*exp(X);
    ST = S(:,end);
end          