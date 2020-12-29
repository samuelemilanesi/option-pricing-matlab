%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price American Put option in the B&S model
% using PDE with theta method on LOG-PRICE transformation
% and SOR algorithm to solve the linear sistem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parameters
% Market parameters
r = 0.01;               % riskfree interest rate 
S0 = 100;               % spot price
% Model parameters
sigma = 0.3;            % standard deviation 
% Contract 
T = 1;                  % maturity
K = S0;                 % strike
payoff = @(x) max(K-S0*exp(x),0);   % payoff function -- Put Option
% Domain Boundaries 
xmax =  (r - sigma^2/2)*T + 6*sigma*sqrt(T);
xmin =  (r - sigma^2/2)*T - 6*sigma*sqrt(T);
upper_boundary_condition = @(t) 0;                                  % v(xmax,t) for EU put
ower_boundary_condition = @(t) K*exp(-r*t)-S0*exp(xmin);            % v(xmin,t) for EU put 
% Numerical schema parameters 
M = 50; dt = T/M;                                                   % time grid
N = 1000; dx = (xmax-xmin)/N; x_grid = linspace(xmin,xmax, N+1)';   % space grid
theta = 0;    % 0.5 = Crank-Nicholson, 0 = Implicit Euler, 1 = Explicit Euler (not uncond stable!)
% SOR params 
tol=1e-4; maxiter=500; omega=1.5;


%% Resolving matrix              
a_d = (1-theta)*(-(r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));    % coeff of v_{i-1,j}
a_m = -1/dt+(1-theta)*(-sigma^2/(dx^2)-r);                   % coeff of v_{i,j} 
a_u = (1-theta)*((r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));     % coeff of v_{i+1,j}
Ahat = spdiags([a_d*ones(N-1,1), a_m*ones(N-1,1), a_u*ones(N-1,1)], -1:1, N-1, N-1);
A = sparse(N+1,N+1); A(2:N,2:N)=Ahat; clear Ahat;
A(1,1)=1; A(2,1)=a_d; A(end,end)=1; A(end-1,end)=a_u;
% Constant term at step j 
b_d = -theta*(-(r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));    % coeff of v_{i-1,j+1}
b_m = -1/dt-theta*(-sigma^2/(dx^2)-r);                    % coeff of v_{i,j+1} 
b_u = -theta*((r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));     % coeff of v_{i+1,j+1}
Bhat = spdiags([b_d*ones(N-1,1), b_m*ones(N-1,1), b_u*ones(N-1,1)], -1:1, N-1, N-1);
B = sparse(N+1,N+1); B(2:N,2:N)=Bhat; clear Bhat;
B(2,1)=b_d; B(end-1,end)=b_u; 


%% Backward in time procedure
v = payoff(x_grid);
for j = (M-1):-1:0
    % compute rhs and adjust with bounday conditions
    rhs = B*v;
    rhs(1) = lower_boundary_condition(T-j*dt);
    rhs(end) = upper_boundary_condition(T-j*dt);
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

fprintf('American Put price under BS model using PDE theta method with theta=%d',theta)
american_put_price= interp1(S0*exp(x_grid), v, S0, 'spline')
%plot(S0*exp(x_grid),v); title('Price'); xlabel('S - spot price');
%hold on
%plot(S0*exp(x_grid),max( K-S0*exp(x_grid),0))
