% Computes FTCS on a convection-diffusion problem in 1/2D (3D in the
% future) 
% Created by: Justin Williams; Last Updated November 15, 2017


function [U, V] = FTCS(Finit, x,y, dx, dt, tf, kappa, U, varargin)

	% TESTING, homework 8b
	[U, V] = FTCSTest(Finit, x,y, dx, dt, tf, kappa, U);
	return;
	
	
	dim = size(Finit, 1);
	y = Finit;
	if(length(varargin)==1)
		if(dim==1)
			y = FTCSNonLinear(Finit, dx, dt, tf, kappa, U);
			return;
		elseif(dim==2)
			y = FTCS2NL(Finit, dx, dt, tf, kappa, U);
			return;
		else
			printc('ERROR: Only 1/2D Nonlinear FTCS implemented!\n');
			return;
		end
	end
	if(dim==1)
		y = FTCS1D(Finit, dx, dt, tf, kappa, U);
	elseif (dim==2)
		y = FTCS2D(Finit,dx,dt,tf,kappa,U);
	elseif (dim == 3)
		fprintf('ERROR: 3D FTCS is not implemented.\n');
	end
end

%%%
% <latex>
% $$
% u_{i,j}^* = u_{i,j}^n - \frac{Co_{x,i,j}}{2} [u_{i+1,j}^n - u_{i-1,j}^n] + \alpha_x [u_{i+1,j}^n - 2 u_{i,j}^n + u_{i-1,j}^n]
% $$
% $$
% v_{i,j}^* = v_{i,j}^n - \frac{Co_x}{2} [v_{i+1,j}^n - v_{i-1,j}^n] + \alpha_x [v_{i+1,j}^n - 2 v_{i,j}^n + v_{i-1,j}^n]
% $$
% $$
% u_{i,j}^{n+1} = u_{i,j}^* - \frac{Co_y}{2} [u_{i,j+1}^* - u_{i,j-1}^*] + \alpha_y [u_{i,j+1}^* - 2 u_{i,j}^* + u_{i,j-1}^*]
% $$
% $$
% v_{i,j}^{n+1} = v_{i,j}^* - \frac{Co_y}{2} [v_{i,j+1}^* - v_{i,j-1}^*] + \alpha_y [v_{i,j+1}^* - 2 v_{i,j}^* + v_{i,j-1}^*]
% $$
% </latex>
%%%
function [U, V] = FTCSTest(Finit, X, Y, dx, dt, tf, kappa, exactFunc) 
	dy = dx(2); dx = dx(1); nx = size(X,1); ny = size(Y,1);
	t = 0;
	U = Finit(:,:,1);
	V = Finit(:,:,1);
	
	[Utop, Vtop] = exactFunc(X(1,:),  Y(1,:),  t);
	[Ubot, Vbot] = exactFunc(X(1,:),  Y(ny,:), t);
	[Ulef, Vlef] = exactFunc(X(:,1),  Y(:,1),  t);
	[Urig, Vrig] = exactFunc(X(:,nx), Y(:,ny), t);
	U(nx,:) = Ubot; V(nx,:) = Vbot;
	U(1,:)  = Utop; V(1,:)  = Vtop;
	U(:,1)  = Ulef; V(:,1)  = Vrig;
	U(:,ny) = Urig; V(:,ny) = Vlef;		
	US = U; VS = V;
		
	while(t<tf+dt)
%		figure(7); surf(X,Y,U); title(t); pause(.5); 
		Co_x = U(2:nx-1,2:ny-1)*dt/dx; Co_y = V(2:nx-1,2:ny-1)*dt/dx;
		alpha_x = kappa*dt/dx^2; alpha_y = kappa*dt/dy^2;
		
		US(2:nx-1,2:ny-1) = U(2:nx-1,2:ny-1) - Co_x/2.*(U(3:nx,2:nx-1)-U(1:nx-2,2:nx-1))+alpha_y*(U(3:nx,2:nx-1) - 2*U(2:nx-1,2:nx-1) + U(1:nx-2,2:nx-1));
		VS(2:nx-1,2:ny-1) = V(2:nx-1,2:ny-1) - Co_x/2.*(V(3:nx,2:nx-1)-V(1:nx-2,2:nx-1))+alpha_y*(V(3:nx,2:nx-1) - 2*V(2:nx-1,2:nx-1) + V(1:nx-2,2:nx-1));
		
		U(2:nx-1,2:ny-1) = US(2:nx-1,2:ny-1) - Co_y/2.*(US(2:nx-1,3:nx)-US(2:nx-1,1:nx-2)) + alpha_x*(US(2:nx-1,3:nx)-2*US(2:nx-1,2:nx-1)+US(2:nx-1,1:nx-2));
		V(2:nx-1,2:ny-1) = VS(2:nx-1,2:ny-1) - Co_y/2.*(VS(2:nx-1,3:nx)-VS(2:nx-1,1:nx-2)) + alpha_x*(VS(2:nx-1,3:nx)-2*VS(2:nx-1,2:nx-1)+VS(2:nx-1,1:nx-2));
		t = t+dt;
	end
end

function y = FTCS1D(Finit, dx, dt, tf, kappa, U) 
	n = size(Finit,2); F1 = Finit; F2=Finit; t = 0;
	alpha = kappa*dt/dx^2; CFL = U*dt/dx; Re = U*dx/kappa;
	if(CFL > 1 || Re > 2)
		fprintf('WARNING: FTCS solution unstable! Use different method.\n')
	end
	
	while(t<tf)
		for i = 2:n-1
			F2(i) = F1(i) - CFL/2 * (F1(i+1) - F1(i-1)); % Convection
			F2(i) = F2(i) + alpha*(F1(i+1) + F1(i-1) - 2*F1(i)); % Diffusion
		end
		F1 = F2; t = t+dt;
	end
	y = F1;
end

function y = FTCS2D(Finit, dx, dt, tf, kappa, U)
% FTCS Performs 2D FTCS scheme on Finit
%	Created 10/26/2017 by Justin Williams
	if(U ~= 0)
		fprintf('ERROR: Convection not implemented in 2D FTCS\n');
		y = Finit;
		return;
	end
	
	dy = dx(2); dx = dx(1);
	
    T2 = zeros(size(Finit));
    T1 = Finit;
    Nx = size(Finit,1);
    Ny = size(Finit,2);
    
    alpha_x = kappa*dt/(dx^2);
    alpha_y = kappa*dt/(dy^2);
    alpha = alpha_x;
    if(alpha_x+alpha_y>1/2) 
        fprintf('WARNING: FTCS in unstable regime\n');
    end
    t = 0;
    while (t < tf) 
        for(i=2:Nx-1)
           for(j=2:Nx-1)
              T2(i,j) = T1(i,j) + alpha*(T1(i+1,j) + T1(i-1,j) + T1(i,j+1) + T1(i,j-1) - 4*T1(i,j));
           end
        end
        T1 = T2;
        t = t+dt;
	end
    fprintf('2D FTCS Finished.\n');
    y = T1;
end

function [y, Re, Co] = FTCSNonLinear(Finit, dx, dt, tf, kappa, U)
	fprintf('Starting Non-Linear FTCS...\n');
	n = size(Finit,2); F1 = Finit; F2=Finit; t = 0;
	alpha = kappa*dt/dx^2; 
	rerror = inf;
	Re = zeros(1,n); Co = zeros(1,n);
	while(t<tf)
		for i = 2:n-1
			CFL = -F1(i)*dt/dx; Rey = -F1(i)*dx/kappa; 
			Re(i) = Rey; Co(i) = CFL;
			F2(i) = F1(i) + CFL/2 * (F1(i+1) - F1(i-1)); % Convection
			F2(i) = F2(i) + alpha*(F1(i+1) + F1(i-1) - 2*F1(i)); % Diffusion
		end
		%BCs
		L=1;
		F2(1) = -2*kappa/(L*exp(-kappa*t/L^2));
		F2(n-1) = -2*kappa*cosh(4/L)/(L*sinh(4/L)+L*exp(-kappa*t/L^2));
		rerror = norm(F1 - F2);
		F1 = F2; t = t+dt;
	end
	fprintf('Estimated error: %e\n', rerror);
	y = F1;
	figure(9); semilogy(0:dx:4, Co, 0:dx:4, Re); legend('Co', 'Re_c')
	fprintf('Finished Non-Linear FTCS\n');
end

function y = FTCS2NL(Finit, dx, dt, tf, kappa, U) 
	exact = U; dy = dx(2); dx = dx(1); 
	nx = size(Finit,1); ny = size(Finit,2);
	y = Finit;
	fprintf('WARNING: 2D nonlinear FTCS not implemented!');
	while(t<tf)
		for i = 1:nx
			for j = 1:ny
				
			end
		end
		t = t+dt;
	end
end