%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put and Options under Kou model
% using PDE with theta method on LOG-PRICE transformation
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
payoff = @(x) max(S0*exp(x)-K,0);   % payoff function -- Call Option

% Domain Boundaries 
xmax =  (r - sigma^2/2)*T + 6*sigma*sqrt(T);
xmin =  (r - sigma^2/2)*T - 6*sigma*sqrt(T);

xmax = log(3); xmin= log(0.2);
upper_boundary_condition = @(t) S0*exp(xmax)-K*exp(-r*t);            % v(xmax,t) for EU call
lower_boundary_condition = @(t) 0;                                   % v(xmin,t) for EU call
% Numerical schema parameters 
M = 75; dt = T/M;                                                   % time grid
N = 2000; dx = (xmax-xmin)/N; x_grid = linspace(xmin,xmax, N+1)';   % space grid
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
a_m = -1/dt+(1-theta)*(-sigma^2/(dx^2)-r-lambdaK);                 % coeff of v_{i,j} 
a_u = (1-theta)*((r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));     % coeff of v_{i+1,j}
A = spdiags([a_d*ones(N+1,1), a_m*ones(N+1,1), a_u*ones(N+1,1)], -1:1, N+1, N+1);
% Constant term at step j 
b_d = -theta*(-(r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));    % coeff of v_{i-1,j+1}
b_m = -1/dt-theta*(-sigma^2/(dx^2)-r-lambdaK);                  % coeff of v_{i,j+1} 
b_u = -theta*((r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));     % coeff of v_{i+1,j+1}
B = spdiags([b_d*ones(N+1,1), b_m*ones(N+1,1), b_u*ones(N+1,1)], -1:1, N+1, N+1);

%% Backward in time procedure
v = payoff(x_grid);
for j = (M-1):-1:1
    % compute the integral
    I=JD_integral(nu,x_grid,v,ynodes,lower_boundary_condition(T-(j+1)*dt), upper_boundary_condition(T-(j+1)*dt) );
    % compute rhs and adjust with bounday conditions
    rhs = B*v - I;
    % rhs(1) = rhs(1)+ lower_boundary_condition(T-j*dt); % lower_bound.. = 0 in call case
    rhs(end) = rhs(end)+upper_boundary_condition(T-j*dt);
    %compute prices at time tj
    v = A\rhs;
end

fprintf('EU Call price under Kou model using PDE theta method with theta=%d and operator splitting for JD processes',theta)
EU_call_price= interp1(S0*exp(x_grid), v, S0, 'spline')
