%% Problem 1
% 1D Explicit FTCS vs. Crank-Nicolsen
clear; clc; clf;
dt = 1.5; t_final=9; dx = 0.05; kappa = 10^-3; alpha = kappa*dt/dx^2;
xMin = -3; xMax = 3; xx = xMin:dx:xMax;

Tinit = 2;      % Initial Temperature step
TinitMin = -1;  % Range of initial temperature step (in x)
TinitMax = 1;   % Range of initial temperature step

% Initial Conditions
N = (xMax - xMin)/ dx + 1;
N_lower = (TinitMin - xMin) ./ dx + 1;
N_upper = (TinitMax - xMin) ./ dx + 1;

T0 = zeros(N,1); T0(N_lower:N_upper) = Tinit;
Tftcs = zeros(N,2); Tftcs(:,1) = T0;
TCN = crankNicolsen(T0, dx, dt, 0, t_final, kappa);

for(i = 0:dt:t_final)
   % FTCS
   for j = 2:N-1
      Tftcs(j,2) = Tftcs(j-1,1)*alpha + (1-2*alpha)*Tftcs(j,1) + alpha*Tftcs(j+1,1); 
   end   
   Tftcs(:,1) = Tftcs(:,2);
end

fx = Tinit/2 * (erf((1-xx)./(2*sqrt(kappa*t_final))) - erf(-(1+xx)/(2*sqrt(kappa*t_final))));
plot(xx, fx, xx, TCN, xx, Tftcs)
legend('Actual', 'Crank-Nicolsen', 'FTCS'); xlabel('x'); ylabel('Temperature')

% Check the errors
Crank_max = norm(reshape(TCN-fx' , [], 1), inf);
FTCS_max = norm(reshape(Tftcs-fx', [], 1), inf);
fprintf('FTCS Max error: %d\n', FTCS_max);
fprintf('Crank-Nicolsen Max Error: %d\n', Crank_max);

%%
clc; clf; clear;
xMin = 0; xMax = 1; yMin = 0; yMax = 1;
kappa = 0.1; 
dx = 0.05; dy = dx; dt = (dx*dy /kappa)/4; 
t0 = 0; tFinal = 80*dt; 
x = xMin:dx:xMax; y = yMin:dy:yMax;

Tinit = zeros(size(x,2), size(y,2));
for(i = 1:size(x,2))
    xv = i*dx-dx;
    for(j=1:size(y,2))
        yv = j*dy-dy;
        Tinit(i,j) = xv*(1-xv^5)*yv*(1-yv); 
    end
end
Tftcs = FTCS(Tinit, dx, dy, dt, t0, tFinal, kappa);

%Exact
trunc = 50;
Tact = zeros(size(Tinit));
Tact0 = Tact;
Ttest=0;
for(n=1:trunc)
    for(m=1:trunc)
        const1 = 120*(-n^4*pi^4*(-1)^n + 12*pi^2*n^2*(-1)^n + 24 + 24*(-1)^(1+n));
        const = -const1*(-2+2*(-1)^m) / (n^7*pi^10*m^3);
        for(ii = 1:size(x,2))
            yv = dx*ii-dx;
            for(jj=1:size(y,2))
                xv = dy*jj-dy;
                Tact(ii,jj) = Tact(ii,jj) + const*sin(n*pi*xv)*sin(m*pi*yv)*exp(-(n^2+m^2)*pi^2*kappa*(tFinal+dt));
				Tact0(ii,jj) = Tact0(ii,jj) + const*sin(n*pi*xv)*sin(m*pi*yv)*exp(-(n^2+m^2)*pi^2*kappa*(0));
            end
        end
       
    end
end
errorSoln = abs(Tact0' - Tinit);
error = abs(Tftcs - Tact');
fprintf('L-inf norm of Solution error: %d\n', norm(reshape(errorSoln, [], 1), inf));
fprintf('L-inf norm of 2D FTCS: %d\n', norm(reshape(error, [], 1), inf));
figure(1)
surf(x,y,Tact')
title('Exact')
figure(2)
surf(x,y,error);
title('error');
hold off
% for wolfram:
% sum of m = 1 to 100 of (sum of n=1 to 100 of (-120*(-n^4*pi^4 * (-1)^n + 12*(pi*n)^2*(-1)^n + 24 + 24*(-1)^(1+n)) * (-2+2*(-1)^m) / (n^7*pi^10*m^3) * sin(n*pi*0.698827)*sin(m*pi*0.698827)))