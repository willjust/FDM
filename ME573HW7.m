%% ME 573, CFD
% Homework 7
% Justin Williams, 11/15/2017
% Solve f_t + U f_x = 0
clc; clf; clear;
%initial numbers and such
dx = 0.05; dt = 0.01; t_f = 0.5; x_min = -5; x_max = 5; U = pi;
x = x_min:dx:x_max; F_ini = erf((1-x)/0.25) - erf(-(x+1)/0.25);

change_x = floor(t_f*U/dx)+1; % number of positions actual shifted
F_act = zeros;
for i = change_x:size(x,2)
	F_act(i) = F_ini(i-change_x+1);
end

%% Solvers
F_ftcs = FTCS(F_ini, dx, dt, t_f, 0, U); 
F_btcs = FTBS(F_ini, dx, dt, t_f, 0, U);
F_cn = crankNicolsen(F_ini, dx, dt, t_f, 0, U);

%% Display the goodies
e_f = abs(F_ftcs - F_act); e_b = abs(F_btcs-F_act); e_c = abs(F_cn-F_act);

fprintf('Max errors\n');
fprintf('FTCS: %d \tFTBS: %d\tCrank-Nicolsen: %d\n',...
	norm(e_f,inf), norm(e_b,inf), norm(e_c,inf));

figure(1); plot(x, F_ftcs, x, F_btcs, x, F_cn, x, F_act); 
legend('FTCS', 'FTBS (upwind)', 'Crank-Nicolsen', 'Actual');
xlabel('Position'); ylabel('Value'); title('Convection')

figure(2); plot(x, e_f, x, e_b, x, e_c);
legend('FTCS', 'FTBS', 'Crank-Nicolsen');
xlabel('Position'); ylabel('Error'); title('Errors');