%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put and Options in the B&S model
% using PDE with theta method on LOG-PRICE transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parameters
% Market parameters
r = 0.01;               % riskfree interest rate 
S0 = 100;            % spot price
% Model parameters
sigma = 0.3;         % standard deviation 
% Contract 
T = 1;                  % maturity
K = S0;                 % strike
payoff = @(x) max(S0*exp(x)-K,0);   % payoff function -- Call Option
% Domain Boundaries 
xmax =  (r - sigma^2/2)*T + 6*sigma*sqrt(T);
xmin =  (r - sigma^2/2)*T - 6*sigma*sqrt(T);
upper_boundary_condition = @(t) S0*exp(xmax)-K*exp(-r*t);            % v(xmax,t) for EU call
lower_boundary_condition = @(t) 0;                                   % v(xmin,t) for EU call
% Numerical schema parameters 
M = 500; dt = T/M;                                                   % time grid
N = 3000; dx = (xmax-xmin)/N; x_grid = linspace(xmin,xmax, N+1)';    % space grid
theta = 0;    % 0.5 = Crank-Nicholson, 0 = Implicit Euler, 1 = Explicit Euler (not uncond stable!)

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
    %compute prices at time tj
    v = A\rhs;
end

fprintf('EU Call price under BS model using PDE theta method with theta=%d',theta)
EU_call_price= interp1(S0*exp(x_grid), v, S0, 'spline')

[C_bls,P_bls]=blsprice(S0,K,r,T,sigma)
