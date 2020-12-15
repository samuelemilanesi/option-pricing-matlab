%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price Barrier Down&Out call option in the B&S model
% using PDE with theta method on LOG-PRICE transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parameters
% Market parameters
r = 0.02;               % riskfree interest rate 
S0 = 218.75;            % spot price
% Model parameters
sigma = 0.2516;         % standard deviation 
% Contract 
T = 1;                  % maturity
K = S0;                 % strike
D = 0.98*S0;            % barrier
payoff = @(x) max(S0*exp(x)-K,0).*(S0*exp(x)>D);   % payoff function -- D&O Call Option
% Domain Boundaries 
xmax =  (r - sigma^2/2)*T + 6*sigma*sqrt(T);
xmin =  log(D/S0);
upper_boundary_condition = @(t) S0*exp(xmax)-K*exp(-r*t);            % v(xmax,t) for D&O call
lower_boundary_condition = @(t) 0;                                   % v(xmin,t) for D&O call
% Numerical schema parameters 
M = 5000; dt = T/M;                                                   % time grid
N = 30000; dx = (xmax-xmin)/N; x_grid = linspace(xmin,xmax, N+1)';    % space grid
theta = 0.5;    % 0.5 = Crank-Nicholson, 0 = Implicit Euler, 1 = Explicit Euler (not uncond stable!)

%% Resolving matrix              
a_d = (1-theta)*(-(r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));    % coeff of v_{i-1,j}
a_m = -1/dt+(1-theta)*(-sigma^2/(dx^2)-r);                   % coeff of v_{i,j} 
a_u = (1-theta)*((r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));     % coeff of v_{i+1,j}
A = spdiags([a_d*ones(N+1,1), a_m*ones(N+1,1), a_u*ones(N+1,1)], -1:1, N+1, N+1);
% Constant term at step j 
b_d = -theta*(-(r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));    % coeff of v_{i-1,j+1}
b_m = -1/dt-theta*(-sigma^2/(dx^2)-r);                    % coeff of v_{i,j+1} 
b_u = -theta*((r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));     % coeff of v_{i+1,j+1}
B = spdiags([b_d*ones(N+1,1), b_m*ones(N+1,1), b_u*ones(N+1,1)], -1:1, N+1, N+1);

%% Backward in time procedure
v = payoff(x_grid);
for j = (M-1):-1:1
    % compute rhs and adjust with bounday conditions
    rhs = B*v;
    % rhs(1) = rhs(1) + lower_boundary_condition(T-j*dt); % lower_bound_cond=0 in CALL case
    rhs(end) = rhs(end)+upper_boundary_condition(T-j*dt);
    %compute prices at time tj
    v = A\rhs;
end

fprintf('Barrier D&O Call price under BS model using PDE theta method with theta=%d',theta)
call_price = interp1(S0*exp(x_grid), v, S0, 'spline')
figure
plot(S0*exp(x_grid),v); title('Price'); xlabel('S - spot price');