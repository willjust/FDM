%% ME 573 HW 8A. 
% Justin Williams, November 22, 2017
U = 0.2; tf = 2.5; kappa = 0.01; dx = 0.05; dt = 0.01; x = -5:dx:5;
tinit = 0.1; c_0 = 2*pi; c_i = pi;


Finit = LinAdvecAnalytical(x,c_0,c_i,kappa,1,U,0.1);

%% Solutions
Fcs = FTCS(Finit, dx, dt, tf, kappa, U);
Fbs = FTBS(Finit, dx, dt, tf, kappa, U);
Fcn = crankNicolsen(Finit, dx, dt, tf, kappa, U);

%% Plots 'n stuff
Fexact = LinAdvecAnalytical(x,c_0,c_i,kappa, 1, U, 2.5);

figure(1); plot(x, Fexact, x, Fcs, 'o'); legend('Exact', 'FTCS')
title('Exact vs FTCS; Linear Advection Diffusion')

figure(2); plot(x, Fexact, x, Fbs, 'o'); legend('Exact', 'Upwind');
title('Exact vs Upwind; Linear Advection Diffusion');

figure(3); plot(x, Fexact, x, Fcn, 'o'); legend('Exact', 'Crank-Nicolsen');
title('Exact vs Crank-Nicolsen; Linear Advection Diffusion');

%% Problem 2
clear; clf; %clf(9);
kappa = 0.02; L = 1; tf = 2*L^2/kappa; dx =0.05; x=0:dx:4;
Finit = -2*kappa/L * cosh(x/L) ./ (sinh(x/L) + 1);
dt = 0.05; 
Fcs = FTCS(Finit, dx, dt, tf, kappa, 1,1); %FTCS(Initial, dx, dt, tf, diffusion, advection, nonlinear)
Fcn = crankNicolsen(Finit, dx, dt, tf, kappa, 1,1);

Fexact = -2*kappa/L * cosh(x/L) ./ (sinh(x/L) + exp(-kappa*tf/L^2));
figure(1); plot(x, Fexact, x, Fcs, 'o'); title('Burgers Equations with FTCS');
figure(2); plot(x, Fexact, x, Fcn, 'o'); title('Burgers Equation with CN');