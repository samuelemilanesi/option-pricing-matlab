%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price Barrier D&O  and option under Merton model
% using PIDE with theta method on LOG-PRICE transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parameters
% Market parameters
r = 0.0367;             % riskfree interest rate 
S0 = 100;               % spot price
% Model parameters
sigma = 0.17801;        % BM vol
delta = sqrt(0.4);      % jump vol
mu    = 0.01;           % jump drift
lambdaK = 0.2;          % n jumps intensity
% Contract parameters
T = 1;                  % maturity
K = S0;                 % strike
D = 90;                 % barrier

% Levy measure
nu=@(y) lambdaK./(delta.*sqrt(2*pi)).*exp(-(y-mu).^2./(2*delta.^2));

% Contract 
T = 1;                  % maturity
K = S0;                 % strike
payoff = @(x) max(S0*exp(x)-K,0).*(S0*exp(x)>D);   % payoff function -- DO Call Option

% Domain Boundaries 
xmax = log(3); xmin= log(D/S0); 

upper_boundary_condition = @(t,x) S0*exp(x)-K*exp(-r*t);           % v(x,t) for x>xmax
lower_boundary_condition = @(t,x) 0;                               % v(x,t) for x<xmin

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
    %compute prices at time tj
    v = A\rhs;
end

fprintf('EU Vanilla prices under Merton model using PDE theta method with theta=%d and operator splitting for JD processes',theta)
call_price= interp1(S0*exp(x_grid), v, S0, 'spline')

