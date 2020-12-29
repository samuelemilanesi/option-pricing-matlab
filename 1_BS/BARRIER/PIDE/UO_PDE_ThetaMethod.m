%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price Barrier Up&Out call option in the B&S model
% using PDE with theta method on LOG-PRICE transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parameters
% Market parameters
r = 0.01;               % riskfree interest rate 
S0 = 100;               % spot price
% Model parameters
sigma = 0.3;               % standard deviation 
% Contract 
T = 1;                  % maturity
K = S0;                 % strike
U = 120;                % barrier
payoff = @(x) max(S0*exp(x)-K,0).*(S0*exp(x)<U);   % payoff function -- U&O Call Option
% Domain Boundaries 
xmin = (r - sigma^2/2)*T - 6*sigma*sqrt(T);
xmax =  log(U/S0);
lower_boundary_condition = @(t) 0;                 % v(xmin,t) for U&O call
upper_boundary_condition = @(t) 0;                 % v(xmax,t) for U&O call
% Numerical schema parameters 
M = 500; dt = T/M;                                                   % time grid
N = 3000; dx = (xmax-xmin)/N; x_grid = linspace(xmin,xmax, N+1)';    % space grid
theta = 0.5;    % 0.5 = Crank-Nicholson, 0 = Implicit Euler, 1 = Explicit Euler (not uncond stable!)

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
    rhs(1) = lower_boundary_condition(T-j*dt); % lower_bound_cond=0 in CALL case
    rhs(end) = upper_boundary_condition(T-j*dt);
    %compute prices at time tj
    v = A\rhs;
end

fprintf('Barrier D&O Call price under BS model using PDE theta method with theta=%d',theta)
call_price = interp1(S0*exp(x_grid), v, S0, 'spline')
figure
plot(S0*exp(x_grid),v); title('Price'); xlabel('S - spot price');