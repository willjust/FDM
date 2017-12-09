function F = FletchersExact(x, y, t) 
%% Analytical solution to steady non-linear diffusion-advection.
% Created: Justin Williams 12/7/2017

%%%%%%%%%%%%%%%%
% Parameters 
%%%%%%%%%%%%%%%%
k = 1/2; nu = 1/5; x0 = 1.0;
a0=0.001*k*exp((1+x0)*k); a1=a0; ymax = 2; a2=0; a3=0;

%%%%%%%%%%%%%%%%
% Function
%%%%%%%%%%%%%%%%
t1 = 2*nu*(a1+a3*y + k*(exp(k*(x-x0)) - exp(-k*(x-x0)))*cos(k*y));
t2 = a0+a1*x+a2*y+a3*x*y+(exp(-k*(x-x0)) + exp(k*(x-x0)))*cos(k*y);
t3 = 2*nu*(a2+a3*x - k*(exp(-k*(x-x0)) + exp(k*(x-x0)))*sin(k*y));
T(1) = t1/t2;
T(2) = t3/t2;
end