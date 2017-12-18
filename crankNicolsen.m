%% Implicit Crank-Nicolsen using Thomas's Alogithm
% currently up to 2D CD problems
% Created by: Justin Williams, November 15, 2017
function [U,V] = crankNicolsen(Finit, X,Y,dx, dt, tf, kappa, U, varargin)
[U,V] = CNTest(Finit, X, Y, dx, dt, tf, kappa, U);
return;
dim = size(Finit, 1);
if(length(varargin) == 1)
	if(dim==1)
		y = crankNicolsenNL(Finit,dx,dt,tf,kappa,U);
	elseif(dim==2)
		y = crankNicolsen2NL(Finit, dx, dt, tf, kappa, U);
	return;
	end
end
if(dim==1)
	y = crankNicolsen1D(Finit, dx, dt, tf, kappa, U);
end
end

function [U,V] = CNTest(Finit, X,Y,dx,dt,tf,kappa,exactFunc)
	dy = dx(2); dx = dx(1); nx = size(X,1); ny = size(Y,1);
	t = 0;
	U = Finit(:,:,1);
	V = Finit(:,:,1);
	b = zeros(nx,1);
	A = zeros(nx);
	
	[Utop, Vtop] = exactFunc(X(1,:),  Y(1,:),  t);
	[Ubot, Vbot] = exactFunc(X(1,:),  Y(ny,:), t);
	[Ulef, Vlef] = exactFunc(X(:,1),  Y(:,1),  t);
	[Urig, Vrig] = exactFunc(X(:,nx), Y(:,ny), t);
	U(nx,:) = Ubot; V(nx,:) = Vbot;
	U(1,:)  = Utop; V(1,:)  = Vtop;
	U(:,1)  = Ulef; V(:,1)  = Vlef;
	U(:,ny) = Urig; V(:,ny) = Vrig;
 		
	while(t<tf+dt)		
		Co_x = U*dt/dx; Co_y = V*dt/dx;
		alpha = kappa*dt/dx^2; 
		rhsc1 = alpha/2+Co_x/4; rhsc2 = 1 - alpha; rhsc3 = alpha/2 - Co_x/4;
		lhsc1 = -rhsc1;	lhsc2 = 1 + alpha; lhsc3 = -rhsc3;
		
		% Loop over y
		for i = 2:ny-1
			
			% Create A (spdiags is doing some weird stuff)
			for j = 2:nx-1
				A(j,j) = lhsc2;
				A(j,j+1) = lhsc3(j+1,i);
				A(j,j-1) = lhsc1(j-1,i);
			end
			A(1,1)= lhsc2; A(1,2) = lhsc3(1,i);
			A(nx,nx) = lhsc2; A(nx, nx-1) = lhsc1(nx,i);
			
			% by
			b(2:nx-1) = rhsc3(3:nx)'.*U(3:nx,i) + rhsc2*U(2:nx-1,i) + rhsc1(1:nx-2)'.*U(1:nx-2,i);
			b(1) = U(1,i)+rhsc1(1)*U(1,i);
			b(nx) = U(nx,i)+rhsc3(nx)*U(1,i);
			temp = thomas(A,b);
			U(2:nx-1,i) = temp(2:nx-1);
			
			%bx
			b(2:nx-1) = rhsc3(3:nx)'.*V(3:nx,i) + rhsc2*V(2:nx-1,i) + rhsc1(1:nx-2)'.*V(1:nx-2,i);
			b(1) = V(1,i) + rhsc3(1)*V(1,i);
			b(nx) = V(nx,i) + rhsc1(nx)*V(nx,i);
			temp = thomas(A,b);
			V(2:nx-1,i) = temp(2:nx-1);
		end
		
		% Loop over y
		rhsc1 = alpha/2 + Co_y/4; rhsc2 = 1-alpha; rhsc3 = alpha/2 - Co_y/4;
		lhsc1 = -rhsc1; lhsc2 = 1+alpha; lhsc3 = -rhsc3;

		for i = 2:ny-1
			% Create A
			for j = 2:nx-1
				A(j,j) = lhsc2;
				A(j,j-1) = lhsc1(j-1,i);
				A(j,j+1) = lhsc3(j+1,i);
			end
			A(1,1)= lhsc2; A(1,2) = lhsc3(1,i);
			A(nx,nx) = lhsc2; A(nx, nx-1) = lhsc1(nx,i);

			b(2:nx-1) = rhsc3(3:nx)'.*U(i,3:nx)' + rhsc2*U(i,2:nx-1)' + rhsc1(1:nx-2)'.*U(i,1:nx-2)';
			b(1) = U(i,1) + rhsc3(1)*Ulef(i);
			b(nx) = U(i,nx) + rhsc1(nx).*Urig(i);
			temp = thomas(A,b);
			U(i, 2:nx-1) = temp(2:nx-1);
			
			b(2:nx-1) = rhsc3(3:nx)'.*V(i,3:nx)' + rhsc2*V(i,2:nx-1)' + rhsc1(1:nx-2)'.*V(i,1:nx-2)';
			b(1) = V(i,1) + rhsc1(1)*V(i,1);
			b(nx) = V(i,nx) + rhsc3(nx).*V(i,nx);
			temp = thomas(A,b);
			V(i, 2:nx-1) = temp(2:nx-1);
		end
		figure(7); surf(X,Y,V); title(t); pause(dt); 

		t = t+dt;
	end
end

%% 1D Implicit Crank-Nicolsen using Thomas' Alogithm
% Solves a diffusion problem using a second order accurate Crank-Nicolsen
% Created by: Justin Williams on October 21, 2017, modified November 15
function y = crankNicolsen1D(Finit, dx, dt, tf, kappa, U) 
    y = Finit; n = size(Finit,2); t = 0;
    
    alpha = kappa*dt/dx^2; CFL = U*dt/dx;

    % Define my diagonal matrix
	c1 = -CFL-2*alpha; c2 = 4*(1+alpha); c3 = CFL-2*alpha;
    e = ones(n-2, 1);
    A = spdiags([c1*e c2*e c3*e], -1:1, n-2,n-2);

    while(t<tf) 
        % Define b
        b = (CFL+2*alpha)*y(1:n-2) + (4-4*alpha)*y(2:n-1) - c3*y(3:n);
		%B.C's
		b(1) = b(1) - c1*y(1);
		b(n-2) = b(n-2)-c3*Finit(n-1);
		
        y(2:n-1) = thomasAlg(A,b);
		t = t+dt;
		plot(y)
    end
end

function y = crankNicolsenNL(Finit, dx, dt, tf, kappa, U)
y = Finit; n  = size(Finit, 2); t = 0;
alpha = kappa*dt*1.35/dx^2;
e = ones(n-2,1);

% For BC, make more generalized %
L = 1;

while(t<tf)
	% Define A
	CFL = U*y(2:n-1)*dt/dx; Re = -y(2:n-1)*dt/kappa; 
	c1 = -CFL-2*alpha; c2 = 4*(1+alpha); c3 = CFL - 2*alpha;
	A = spdiags([c1' c2*e c3'], -1:1, n-2,n-2);
	
	% Define b
	b = (CFL+2*alpha).*y(1:n-2) +(4-4*alpha)*y(2:n-1) - c3.*y(3:n);
	
	%BC's
	y(1) = -2*kappa/(L*exp(-kappa*t/L^2));
	b(1) = b(1) + (-CFL(1) + 2*alpha)*y(1);
	y(n) = (-2*kappa*cosh(4/L)/(L*sinh(4/L)+L*exp(-kappa*t/L^2)));
	b(n-2) = b(n-2) + (2*alpha - CFL(n-2))*y(n);
	
	%solve
	y(2:n-1) = thomasAlg(A,b);
	t = t+dt;
end
figure(8); semilogy(dx:dx:4-dx, -CFL, dx:dx:4-dx, Re);
legend('CFL', 'Re');
end

function y = crankNicolsen2NL(Finit, dx, dt, tf, kappa, exact)
y = Finit; nx  = size(Finit, 1); ny = size(Finit, 2); t = 0;
alpha = kappa*dt/dx^2;
e = ones(n-2,1);

while(t<tf)
	for i = 1:Nx
	% Define A
	Co = y(i, :)*dt/dx;
	c1 = -Co - 2*alpha; c2 = 4*(1+alpha); c3 = Co - 2*alpha;
	A = spdiags([c1' c2*e, c2'], -1:1, n-2, n-2);
	disp(full(A));
	
	% Define b

	
	% Solve
	
	end
	for j = 1:Ny
	% Define A
	
	% Define b
	
	% Solve
	
	end
	t = t+dt;
end
end