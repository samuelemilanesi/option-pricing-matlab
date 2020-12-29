%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price American Put option under KOU model
% using PDE with theta method on LOG-PRICE transformation
% and SOR algorithm to solve the linear sistem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parameters
% Market parameters
r = 0.01;               % riskfree interest rate 
S0 = 100;            % spot price
% Model parameters
sigma = 0.3;         % standard deviation 
p = 0.5;             % pos jump prob
lambdaK = 2;         % jump time intensity
lambdap = 15;
lambdam = 18;

% Levy measure
nu=@(y) lambdaK*(p*lambdap*exp(-lambdap*y).*(y>0)+(1-p)*lambdam*exp(-lambdam*abs(y)).*(y<0)); 

% Contract 
T = 1;                  % maturity
K = S0;                 % strike
payoff = @(x) max(K-S0*exp(x),0);   % payoff function -- Put Option
% Domain Boundaries 
xmax = log(3); xmin= log(0.2);

upper_boundary_condition = @(t,x) 0;                               % v(x,t) for x>xmax
lower_boundary_condition = @(t,x) K-S0*exp(x);            % v(x,t) for x<xmin
% Numerical schema parameters 
M = 50; dt = T/M;                                                   % time grid
N = 1000; dx = (xmax-xmin)/N; x_grid = linspace(xmin,xmax, N+1)';   % space grid
theta = 0;    % 0.5 = Crank-Nicholson, 0 = Implicit Euler, 1 = Explicit Euler (not uncond stable!)


%% Computing alpha, and truncating the integral
% Truncate the integral
tol=1e-14;
ymin=-0.1;
while nu(ymin)>tol
    ymin=ymin-0.5;
end
ymax=0.1;
while nu(ymax)>tol
    ymax=ymax+0.5;
end
ynodes=linspace(ymin,ymax,2*N);
%figure
%plot(ynodes,nu(ynodes));
alpha=trapz(ynodes, (exp(ynodes)-1).*nu(ynodes))
%lambdaNum=trapz(ynodes, nu(ynodes)) % rmk: this is an approx of the real lambdaK passed as param

%% Resolving matrix              
a_d = (1-theta)*(-(r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));    % coeff of v_{i-1,j}
a_m = -1/dt+(1-theta)*(-sigma^2/(dx^2)-r-lambdaK);                   % coeff of v_{i,j} 
a_u = (1-theta)*((r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));     % coeff of v_{i+1,j}
Ahat = spdiags([a_d*ones(N-1,1), a_m*ones(N-1,1), a_u*ones(N-1,1)], -1:1, N-1, N-1);
A = sparse(N+1,N+1); A(2:N,2:N)=Ahat; clear Ahat;
A(1,1)=1; A(2,1)=a_d; A(end,end)=1; A(end-1,end)=a_u;
% Constant term at step j 
b_d = -theta*(-(r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));    % coeff of v_{i-1,j+1}
b_m = -1/dt-theta*(-sigma^2/(dx^2)-r-lambdaK);                    % coeff of v_{i,j+1} 
b_u = -theta*((r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));     % coeff of v_{i+1,j+1}
Bhat = spdiags([b_d*ones(N-1,1), b_m*ones(N-1,1), b_u*ones(N-1,1)], -1:1, N-1, N-1);
B = sparse(N+1,N+1); B(2:N,2:N)=Bhat; clear Bhat;
B(2,1)=b_d; B(end-1,end)=b_u; 

% SOR params 
tol=1e-6; maxiter=500; omega=1.5;

%% Backward in time procedure
v = payoff(x_grid);
for j = (M-1):-1:0
    % compute boundary conditions 
    lbc_old = @(x) lower_boundary_condition(T-(j+1)*dt, x);  
    ubc_old = @(x) upper_boundary_condition(T-(j+1)*dt, x);
    lbc = lower_boundary_condition(T-j*dt,xmin);  
    ubc = upper_boundary_condition(T-j*dt,xmax);
    % compute the integral
    I = JD_integral(nu,x_grid,v,ynodes,lbc_old, ubc_old);
    % compute rhs and adjust boundary conditions
    rhs = B*v - I;
    rhs(1) = lbc;
    rhs(end) = ubc;    

    %compute prices at time tj -- exploit SOR to solve the linear system
    % guess solution: v 
    for iter = 1:maxiter 
        v_old = v;          % next iteration guess is the previous result 
        for i=1:length(v)
            if i==1
                iter_rhs = rhs(i)-A(i,i+1)*v_old(i+1);
            elseif i==length(v)
                iter_rhs = rhs(i)-A(i,i-1)*v(i-1);
            else
                iter_rhs = rhs(i)-A(i,i+1)*v_old(i+1)-A(i,i-1)*v(i-1);
            end
            v(i)=max( (1-omega)*v_old(i)+omega*iter_rhs/A(i,i), K-S0*exp(x_grid(i)) );
        end
        err=norm(v-v_old,'Inf');
        if err<tol
            break
        end
    end
    [j, iter, err]
end

fprintf('American Put price under Kou model using PDE theta method with theta=%d',theta)
Kou_american_put_price= interp1(S0*exp(x_grid), v, S0, 'spline')
%plot(S0*exp(x_grid),v); title('Price'); xlabel('S - spot price');
%hold on
%plot(S0*exp(x_grid),max( K-S0*exp(x_grid),0))
