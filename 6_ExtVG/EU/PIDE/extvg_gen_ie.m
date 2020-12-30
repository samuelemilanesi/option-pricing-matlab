clear
close
% PRICE AN EUROPEAN CALL OPTION
% FD SCHEME - IMPLICIT EULER
% LOG-PRICE TRANSFORM
% GENERAL LEVY - EXTENDED VG MODEL
clear all
close all
%% Input
T=1; K=100; S0=100; r=0.01; 
%%%%%%%%%%%%%%%%%%%
% change 1--> VG
theta=0.05; sigma=0.6; k=0.2;  sigmaGBM=0.4;
A=theta/(sigma^2); B=sqrt(theta^2+2*sigma^2/k)/(sigma^2);
nu=@(y) 1./(k*abs(y)).*exp(A*y-B*abs(y));
%%%%%%%%%%%%%%%%%%%
M=50; N=1000; epsilon=0.;
if epsilon==0
    sigma=0;
else
    %%%%%%%%%%%%%%%%%%%
    %  change 2: define sigma from truncation, and define the truncated nu
    ynodes=linspace(-epsilon,epsilon,2*N);
    sigma=sqrt(trapz(ynodes,ynodes.^2.*nu(ynodes)))
    nu=@(y) nu(y).*(abs(y)>epsilon);
    %%%%%%%%%%%%%%%%%%%
end
sigma=sqrt( sigmaGBM^2+sigma^2);
%% Grids
dt=T/M; 
Smin=0.2*S0; Smax=3*S0; xmin=log(Smin/S0); xmax=log(Smax/S0);
dx=(xmax-xmin)/N;
x=linspace(xmin, xmax, N+1)';

%% Truncate the integral
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
figure
plot(ynodes,nu(ynodes));

%% Construct the matrix
A=-(r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2);
B=-1/dt-sigma^2/(dx^2)-(r);
C=(r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2);
% 1nd way
Mhat=spdiags( [A*ones(N-1,1) B*ones(N-1,1) C*ones(N-1,1)],[-1 0 1],N-1,N-1);
% 2nd way
Mat=sparse(N+1,N+1); Mat(2:end-1,2:end-1)=Mhat; clear Mhat
Mat(1,1)=1; Mat(2,1)=A; Mat(end,end)=1; Mat(end-1,end)=C;

%% Backward-in-time procedure
V=max( S0*exp(x)-K,0); 
for j=M-1:-1:0
    I=integral(nu,x,V,ynodes,S0,K*exp(-r*(T-(j+1)*dt)));
    rhs=-1/dt*V-I; rhs(1)=0; rhs(end)=S0*exp(xmax)-K*exp(-r*(T-j*dt));
    V=Mat\rhs;
end
figure
plot(S0*exp(x),V); title('Price'); xlabel('S - spot price');
Price=interp1( S0*exp(x),V, S0,'spline')
Price_CM=FFT_CM_Call5

function I=integral(nu,x,V,ynodes,S0,discK)
% trapezoidal formula
dy=ynodes(2)-ynodes(1);
w=ones(size(ynodes)); w(1)=1/2; w(end)=1/2;
% construct the first order derivative
dx=x(2)-x(1);
dV=(V(3:end)-V(1:end-2))/(2*dx);
% compute the integral
I=zeros(size(x));
for i=2:length(I)-1
    I(i)=sum( w.*(valueV(x(i)+ynodes,x,V,S0,discK)-V(i)-(exp(ynodes)-1)*dV(i-1) )...
             .*nu(ynodes) )*dy;
end
end

function f=valueV(y,x,V,S0,discK)
f=zeros(size(y));
index=find( (y>=x(1)) .* (y<=x(end)) );
f(index)=interp1(x,V,y(index));
index=find( (y>x(end)) );
f(index)=S0*exp(y(index))-discK;
end    
    

    
    
    
    
    
    
    
    
    
    


