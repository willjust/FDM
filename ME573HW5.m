% ME573HW5.m
%% Question 1.
% FTCS Heat Diffusion vs analytical at t=1000
clear; clc; clf;
% Parameters
dx = [0.5, 0.25, 0.05]; % Range of spatial steps
dt = 0.01;      % Timestep
kappa = 10^-3;  % Diffusivity
t_final = 10^2; % Time final
xMin = -3;      % Domain beginning
xMax = 3;       % Domain end
Tinit = 2;      % Initial Temperature step
TinitMin = -1;  % Range of initial temperature step (in x)
TinitMax = 1;   % Range of initial temperature step

% Initial Conditions
N = (xMax - xMin) ./ dx + 1;
alpha = kappa.*dt./dx.^2;
N_lower = (TinitMin - xMin) ./ dx + 1;
N_upper = (TinitMax - xMin) ./ dx + 1;

T5 = zeros(N(1), 2);
T25 = zeros(N(2), 2);
T05 = zeros(N(3), 2);

T5(N_lower(1):N_upper(1), 1) = Tinit;
T25(N_lower(2):N_upper(2), 1) = Tinit;
T05(N_lower(3):N_upper(3), 1) = Tinit;

for i = 1:t_final/dt 
    % Apply Explicit FTCS to all timesteps 
   for j = 2:N(1)-1
      T5(j,2) = T5(j-1,1)*alpha(1) + (1-2*alpha(1))*T5(j,1) + alpha(1)*T5(j+1,1); 
   end
    
   for j = 2:N(2)-1
      T25(j,2) = T25(j-1,1)*alpha(2) + (1-2*alpha(2))*T25(j,1) + alpha(2)*T25(j+1,1); 
   end
    
   for j = 2:N(3)-1
      T05(j,2) = T05(j-1,1)*alpha(3) + (1-2*alpha(3))*T05(j,1) + alpha(3)*T05(j+1,1); 
   end
    
   % Move T[i,2] to T[i,1] 
   T5(:,1) = T5(:,2);
   T25(:,1) = T25(:,2);
   T05(:,1) = T05(:,2);
   
end
xx1 = xMin:dx(1):xMax;
xx2 = xMin:dx(2):xMax;
xx3 = xMin:dx(3):xMax;

xx = xMin:0.01:xMax;
fx = Tinit/2 * (erf((1-xx)./(2*sqrt(kappa*t_final))) - erf(-(1+xx)/(2*sqrt(kappa*t_final))));
plot(xx1, T5, '*', xx2, T25, 'x', xx3, T05, 'o',xx,fx)
legend('dx = 0.05', 'dx = 0.25', 'dx = 0.5', 'Exact')
xlabel('x')
ylabel('Temperature')
title('Problem 1: Various spatial discretizations vs actual')
% TODO add analytical solution
%% Problem 2
% BTCS vs FTCS vs Analytical
clear; clc; clf;
dt = 2;
t_final=10^2;
dx = 0.05;
kappa = 10^-3;
alpha = kappa*dt/dx^2;
xMin = -3;
xMax = 3;
xx = xMin:dx:xMax;
Tforward = zeros(121, 2);
Tbackward = zeros(121,2);

Tinit = 2;      % Initial Temperature step
TinitMin = -1;  % Range of initial temperature step (in x)
TinitMax = 1;   % Range of initial temperature step

% Initial Conditions
N = (xMax - xMin)/ dx + 1;
N_lower = (TinitMin - xMin) ./ dx + 1;
N_upper = (TinitMax - xMin) ./ dx + 1;

Tforward(N_lower:N_upper,1) = Tinit;
Tbackward(N_lower:N_upper,1) = Tinit;

for(i = 0:dt:t_final)
   % FTCS
   for j = 2:N-1
      Tforward(j,2) = Tforward(j-1,1)*alpha + (1-2*alpha)*Tforward(j,1) + alpha*Tforward(j+1,1); 
   end   
   Tforward(:,1) = Tforward(:,2);
end

e = ones(N-2,1);
A = spdiags([-alpha*e (1+2*alpha)*e -alpha*e],-1:1,N-2,N-2);
At = full(A);
for(i=0:dt:t_final)
    % b_1, b_N don't need updating as boundaries are assumed fixed at 0
   Tbackward(2:N-1,2) = thomasAlg(A,Tbackward(2:N-1,1));
   Tbackward(2:N-1,1)=Tbackward(2:N-1,2);
end
xx = xMin:0.05:xMax;
fx = Tinit/2 * (erf((1-xx)./(2*sqrt(kappa*t_final))) - erf(-(1+xx)/(2*sqrt(kappa*t_final))));
plot(xx, Tforward,'o', xx,Tbackward, '*', xx, fx)
legend('Actual', 'BTCS *', 'FTCS')
title('Stability of BTCS vs FTCS')
xlabel('x')
ylabel('Temperature')
%axis([xMin xMax 0 2])