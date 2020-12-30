%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price American Put option under VG model
% using PDE with theta method on LOG-PRICE transformation
% and SOR algorithm to solve the linear sistem.
% General Lévy operator splitting. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parameters
% Market parameters
r = 0.01;            % riskfree interest rate 
S0 = 100;            % spot price
% Model parameters
sigma = 0.6;
theta = 0.05;
kVG   = 0.2;
A=theta/(sigma^2); B=sqrt(theta^2+2*sigma^2/kVG)/(sigma^2);
% Contract parameters
T = 1;                  % maturity
K = S0;  % strike

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'theta',theta,'kVG',kVG);

% Levy measure
nu=@(y) 1./(kVG*abs(y)).*exp(A*y-B*abs(y));
% Contract 
T = 1;                  % maturity
K = S0;                 % strike
payoff = @(x) max(K-S0*exp(x),0);   % payoff function -- Put Option
% Domain Boundaries 
xmax = log(3); xmin= log(0.2);

upper_boundary_condition = @(t,x) 0;                       % v(x,t) for x>xmax
lower_boundary_condition = @(t,x) K-S0*exp(x);            % v(x,t) for x<xmin
% Numerical schema parameters 
M = 50; dt = T/M;                                                   % time grid
N = 1000; dx = (xmax-xmin)/N; x_grid = linspace(xmin,xmax, N+1)';   % space grid
theta = 0;    % 0.5 = Crank-Nicholson, 0 = Implicit Euler, 1 = Explicit Euler (not uncond stable!)

%% Truncating small jumps and evaluating sigma of characteristic triplet according to A-R and truncated Lévy measure

%{
  rmk: the truncation of small jumps is performed because of
  the numerical schema we use is for parabolic PDE and if sigma=0
  this is no more true. Change the schema one could avoid small jumps truncation
%}
epsilon = 0.1;
if epsilon==0
    sigma = 0;
else
    ynodes=linspace(-epsilon,epsilon,2*N);
    sigma=sqrt(trapz(ynodes,ynodes.^2.*nu(ynodes)));
    nu=@(y) nu(y).*(abs(y)>epsilon);
end

%% Truncating the integral
tol=1e-14;
ymin=-0.1-epsilon;
while nu(ymin)>tol
    ymin=ymin-0.5;
end
ymax=0.1+epsilon;
while nu(ymax)>tol
    ymax=ymax+0.5;
end
ynodes=linspace(ymin,ymax,2*N);

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
    I = GenLevy_integral(nu,x_grid,v,ynodes,lbc_old, ubc_old);
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

fprintf('American Put price under VG model using PIDE theta method with theta=%d and operator splitting for general Lévy but truncating small jumps',theta)
VG_american_put_price= interp1(S0*exp(x_grid), v, S0, 'spline')
plot(S0*exp(x_grid),v); title('Price'); xlabel('S - spot price');
hold on
plot(S0*exp(x_grid),max( K-S0*exp(x_grid),0))
