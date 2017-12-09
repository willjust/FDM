%% ME 573 HW 6 
% Justin Williams, November 5, 2017
clear; clc; clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%						Setup Parameters								  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin = 0; xmax = 1; ymin = 0; ymax = 1;
dx = 0.01; dy = dx;
nx = (xmax-xmin)/dx + 1; ny = (ymax-ymin)/dy + 1;
omega = 1.8; omega_opt = 2/(1+ sin(pi/(nx+1)));
f_0 = zeros(nx,ny);
f_actual = f_0;
for i = 2:nx-1
	x = i*dx - dx;
	for j = 2:ny-1
		y = j*dy - dy;
		f_0(i,j) = -2*x*(1-x) - 2*y*(1-y);
		f_actual(i,j) = x*(1-x)*y*(1-y);
	end
end

iteration = 800;

f_jacobi = zeros(size(f_0));
f_tmp = f_jacobi; f_gs = f_tmp; f_sor2 = f_tmp; f_sor1 = f_tmp; 
f_sur1 = f_tmp; f_sur2 = f_tmp; f_sor4 = f_tmp; f_sor3 = f_tmp; 
error = zeros(iteration, 5); 
for t=1:iteration
	for i=2:nx-1
		x = i*dx-dx;
		for j = 2:ny-1
			y = j*dy - dy;
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%					Jocobi Method							  %
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			f_tmp(i,j) = 1/4*(f_jacobi(i+1,j) + f_jacobi(i-1,j)...
				+ f_jacobi(i,j+1) + f_jacobi(i,j-1))...
				-1/4*dx^2*(-2*x*(1-x)-2*y*(1-y));
			error(t,1) = norm(reshape(f_actual-f_tmp,1,[]),inf);
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%					Gauss-Seidel		   					  %	
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			f_gs(i,j) = 1/4*(f_gs(i-1,j) + f_gs(i+1,j) + f_gs(i,j+1)...
				+ f_gs(i,j-1)) - 1/4* (-2*x*(1-x) - 2*y*(1-y))*dx^2;
			error(t,2) = norm(reshape(f_actual-f_gs,1,[]),inf);
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%					SOR Method								  %
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			% Keep in mind, sor1 has a mix of this time level and next
			f_sor2(i,j) = 1/4*(f_sor1(i+1,j) + f_sor1(i-1,j)...
				+ f_sor1(i,j+1) + f_sor1(i,j-1))...
				-1/4*dx^2*(-2*x*(1-x)-2*y*(1-y));
			f_sor1(i,j) = f_sor1(i,j) + omega*(f_sor2(i,j)-f_sor1(i,j));
			error(t,3) = norm(reshape(f_actual-f_sor1,1,[]),inf);
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%					SOR-Optimal Method						  %
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			% Keep in mind, sor3 has a mix of this time level and next
			f_sor4(i,j) = 1/4*(f_sor3(i+1,j) + f_sor3(i-1,j)...
				+ f_sor3(i,j+1) + f_sor3(i,j-1))...
				-1/4*dx^2*(-2*x*(1-x)-2*y*(1-y));
			f_sor3(i,j) = f_sor3(i,j) + omega*(f_sor4(i,j)-f_sor3(i,j));
			error(t,5) = norm(reshape(f_actual-f_sor3,1,[]),inf);
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%					SUR Method								  %
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			% For fun, what happens when \omega < 1; here it is 0.5
			f_sur2(i,j) = 1/4*(f_sur1(i+1,j) + f_sur1(i-1,j)...
				+ f_sur1(i,j+1) + f_sur1(i,j-1))...
				-1/4*dx^2*(-2*x*(1-x)-2*y*(1-y));
			f_sur1(i,j) = f_sur1(i,j) + 0.5*(f_sur2(i,j)-f_sur1(i,j));
			error(t,4) = norm(reshape(f_actual-f_sur1,1,[]),inf);
		end
	end
	f_jacobi = f_tmp;
end
%% Data Processing
e_j = f_actual - f_jacobi; e_gs = f_actual-f_gs; e_sor = f_actual-f_sor1;
e_sur = f_actual-f_sur1; e_soropt = f_actual - f_sor3;

max_j = norm(reshape(e_j,[],1),2); max_gs = norm(reshape(e_gs,[],1),2);
max_sor = norm(reshape(e_sor,[],1),2); max_sur = norm(reshape(e_sur,[],1),2);
max_soropt = norm(reshape(e_soropt, [],1),2);

fprintf('Jacobi 2-norm: %d\n', max_j);
fprintf('Gauss-Seidel 2-norm: %d\n', max_gs);
fprintf('SOR 2-norm (omega=1.8): %d\n', max_sor);
fprintf('SOR 2-norm (omega=%d): %d\n', omega_opt, max_sor);
fprintf('SUR 2-norm (omega=0.5): %d\n', max_sur);

semilogy(error); title('Error of Fixed Point Iterative Methods, dx = 0.01');
xlabel('Iteration'); ylabel('Error'); 
legend('Jacobi', 'Gauss-Seidel', 'SOR', 'SUR', 'SOR Optimal w');