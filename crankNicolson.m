function [U, V] = crankNicolson(U,V,dx,dy, dt, tf, kappa)
t = 0;
[nx, ny] = size(U);
%% Remove X/Y, for visualization
for j = 1:ny % y first to maintain spatial locality of matlab (col-major order)
	for i = 1:nx
		X(i,j) = i*dx - dx; % not using meshgrid as it isn't too intuitive
		Y(i,j) = j*dy - dy; % this orders X,Y as I want them.
	end
end


while t < tf
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%					Convection-Diffusion (Crank-Nicolson)				       %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Constant values Values
	Cox = U*dt/dx; Coy = V*dt/dy; 
	alphax = kappa*dt/dx^2; alphay = kappa*dt/dy^2;

	% Form A_x^n
	A = zeros(nx,nx);
	for i = 2:ny-1		
		A(i, i+1) = Cox(i,i+1)/4 - alphax/2;
		A(i,i)	   = 1 + alphax;
		A(i, i-1) = -Cox(i, i-1)/4 - alphax/2;
	end
	A(1,1) = 1+alphax;   A(1,2) =     Cox(1,2)/4 - alphax/2; 
	A(nx,nx) = 1+alphax; A(nx, nx-1) = -Cox(nx,nx-1)/4 - alphax/2; 

	
	%% REMOVE
	const = 1;
	
	
	%loop over x values
	for i = 2:ny-1
		bx = zeros(nx,1);
		by = zeros(nx, 1);
		bx(1) = U(1,i) + const*(Cox(1,i)/4+alphax/2)*U(1,i);
		by(1) = V(1,i) + const*(Cox(1,i)/4+alphax/2)*V(1,i);
		bx(2:nx-1) = (Cox(1:nx-2,i)/4 + alphax/2).*U(1:nx-2,i)+ (1- alphax)*U(2:nx-1,i)+ (-Cox(3:nx,i)/4 + alphax/2).*U(3:nx,i);
		by(2:nx-1) = (Cox(1:nx-2,i)/4 + alphax/2).*V(1:nx-2,i)+ (1- alphax)*V(2:nx-1,i)+ (-Cox(3:nx,i)/4 + alphax/2).*V(3:nx,i);
		bx(nx) = U(nx,i) + const*(-Cox(nx,i)/4+alphax/2)*U(nx,i);
		by(nx) = V(nx,i) + const*(-Cox(nx,i)/4+alphax/2)*V(nx,i);

		tempx = thomas(A,bx);
		tempy = thomas(A,by);
		
		U(2:nx-1,i) = tempx(2:nx-1);
		V(2:nx-1,i) = tempy(2:nx-1);
	end

	% Form A_y^*
	A = zeros(ny,ny);
	for i = 2:ny-1		
		A(i, i+1) = Coy(i,i+1)/4 - alphay/2;
		A(i,i)	   = 1 + alphay;
		A(i, i-1) = -Coy(i, i-1)/4 - alphay/2;
	end
	A(1,1) = 1+alphay;   A(1,2) =     Coy(1,2)/4 - alphay/2; 
	A(nx,nx) = 1+alphay; A(nx, nx-1) = -Coy(nx,nx-1)/4 - alphay/2; 

	%loop over y values
	for i = 2:nx-1
		bx = zeros(ny,1);
		by = zeros(ny, 1);

		bx(1) = U(i,1) + const*(Coy(i,1)/4+alphay/2)*U(i,1);
		by(1) = V(i,1) + const*(Coy(i,1)/4+alphay/2)*V(i,1);
		bx(2:ny-1) = (Coy(i,1:ny-2)/4 + alphay/2).*U(i,1:ny-2)+ (1-alphay)*U(i,2:ny-1)+ (-Coy(i,3:ny)/4 + alphay/2).*U(i,3:ny);
		by(2:ny-1) = (Coy(i,1:nx-2)/4 + alphay/2).*V(i,1:nx-2)+ (1-alphay)*V(i,2:ny-1)+ (-Coy(i,3:nx)/4 + alphay/2).*V(i,3:nx);
		bx(ny) = U(i,ny) + const*(-Coy(i,nx)/4+alphay/2)*U(i,nx);
		by(ny) = V(i,ny) + const*(-Coy(i,nx)/4+alphay/2)*V(i,nx);

		tempx = thomas(A,bx);
		tempy = thomas(A,by);
		U(i,2:nx-1) = tempx(2:ny-1);
		V(i,2:nx-1) = tempy(2:ny-1);	
	end
	%figure(7); surf(X,Y,V); pause(dt); 
	t=t+dt;
end

end