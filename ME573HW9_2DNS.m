% 2D Navier Stokes Solver
clear;

% Parameters
rho = 1.0; nu = 0.01;
kappa = nu/rho; %nu/rho; laplacian term constant
t0 = 0; dt = 0.1; tfinal = dt; % time parameters
Lx = 1.0; Ly = Lx; dx = 0.02; dy = dx; % Domain boundaries
Vlid = 1.0; sliptol = 1E-3; nslip = 50; % Boundary conditions

% Pressure
load('converged_phi.mat')
% eq = @(x,y) x.*y + 2*x + 2*y;

% phipx = @(x,y) y+2;
% phipy = @(x,y) x+2;

% Domain Creation
xp = 0:dx:Lx; yp = 0:dx:Ly;
nx = size(xp,2); ny = size(yp,2);
U = zeros(nx, ny); V = zeros(nx,ny); X = U; Y = U;

for j = 1:ny % y first to maintain spatial locality of matlab (col-major order)
	for i = 1:nx
		X(i,j) = i*dx - dx; % not using meshgrid as it isn't too intuitive
		Y(i,j) = j*dy - dy; % this orders X,Y as I want them.
	end
end

% Setup initial conditions
for i = 1:nx
	U(i, ny) = Vlid;
end
% surf(X,Y,U); xlabel('X'); ylabel('Y'); % BC Verified
 t = 0;
while(t<tfinal)
maxslip = 5; iteration = 0;
USS = U; VSS = V;
while(maxslip>sliptol)
	maxslip = 0;
	[US,VS] = crankNicolson(USS,VSS,dx,dy,dt,dt,kappa);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%								Projection step                                %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% phi = eq(X,Y);
	PPX = zeros(nx,ny); %Phi prime x
	PPY = zeros(nx,ny); %Phi prime y

	% internal nodes
	for j = 2:ny-1
		for i = 2:nx-1
			PPX(i,j) = (phi(i+1,j) - phi(i-1,j))/2/dx;
			PPY(i,j) = (phi(i,j+1) - phi(i,j-1))/2/dy;
		end
	end
	% Boundary nodes

	for i=3:nx-2
		% j=1 (Bottom nodes)
		PPX(i,1) = (phi(i+1,1)-phi(i-1,1))/2/dx;
		PPY(i,1) = 0;
		% j=ny (Top nodes)
		PPX(i,ny) = (phi(i+1,ny)-phi(i-1,ny))/2/dx;
		PPY(i,ny) = 0;
	end
	PPX(2,1) = (-3*phi(2,1) + 4*phi(3,1) - phi(4,1))/2/dx; % Bottom Left
	PPX(2,ny) = (-3*phi(2,ny)+4*phi(3,ny) - phi(4,ny))/2/dx; % Bottom Right
	PPX(nx-1,1) = (3*phi(nx-1,1)-4*phi(nx-2,1) + phi(nx-3,1))/2/dx; % Top Left
	PPX(nx-1,ny) = (3*phi(nx-1,ny)-4*phi(nx-2,ny) + phi(nx-3,ny))/2/dx; % Top Right

	for j=2:ny-1
		% i = 1 (Left nodes)
		PPX(1,j) = 0;
		PPY(1,j) = (phi(1,j+1)-phi(i,j-1))/2/dy;
		% i = nx (Right nodes)
		PPX(nx,j) = 0;
		PPY(nx,j) = (phi(nx,j+1)-phi(i,j-1))/2/dy;
	end
	PPY(1,2) = (-3*phi(1,2) + 4*phi(1,3) - phi(1,4))/2/dx; % Bottom Left
	PPY(nx,2) = (-3*phi(nx,2)+4*phi(nx,3) - phi(nx,4))/2/dx; % Bottom Right
	PPY(1,ny-1) = (3*phi(1,ny-1)-4*phi(1,ny-2) + phi(1,ny-3))/2/dx; % Top Left
	PPY(nx,ny-1) = (3*phi(nx,ny-1)-4*phi(nx,ny-2) + phi(nx,ny-3))/2/dx; % Top Right

	%Pactx = phipx(X,Y); actual d(Phi)/dx
	%Pacty = phipy(X,Y); actual d(Phi)/dy

	% Actually do the projection
	% Internal Update
	US(2:nx-1,2:ny-1) = US(2:nx-1,2:ny-1)-dt/rho*PPX(2:nx-1,2:ny-1);
	VS(2:nx-1,2:ny-1) = VS(2:nx-1,2:ny-1)-dt/rho*PPY(2:nx-1,2:ny-1);
	% Left Edge
	US(1,:) = 0 - dt/rho*PPX(1,:);
	VS(1,:) = 0 - dt/rho*PPY(1,:);
	% Right Edge
	US(nx,:) = 0-dt/rho*PPX(nx,:);
	VS(nx,:) = 0-dt/rho*PPY(nx,:);
	% Bottom Edge
	US(:,1) = 0-dt/rho*PPX(:,1);
	VS(:,1) = 0-dt/rho*PPX(:,1);
	% Top Edge (Vlid)
	US(:,ny) = Vlid - dt/rho*PPX(:,ny);
	VS(:,ny) = 0 - dt/rho*PPX(:,ny);
	
	% Calculate slip
	% Left Edge
	slip = VS(1,:);
	VSS(1,:) = VSS(1,:) - slip;
	slip=max(abs(slip));
	if(slip>maxslip) 
		maxslip = slip;
	end
	
	% Right Edge
	slip = VS(nx,:);
	VSS(nx,:) = VSS(nx,:) - slip;
	slip=max(abs(slip));
	if(slip>maxslip)
		maxslip = slip;
	end
	
	% Bottom Edge
	slip = US(:,1);
	USS(:,1) = USS(:,1) - slip;
	slip=max(abs(slip));
	if(slip>maxslip)
		maxslip=slip;
	end
	
	% Top Edge
	slip = US(:,ny)-Vlid;
	USS(:,ny) = USS(:,ny) - slip;
	slip=max(abs(slip));
	if(slip>maxslip)
		maxslip=slip;
	end
	iteration = iteration + 1;
	figure(1); surf(X,Y,USS); title(iteration); pause(dt);
	figure(2); surf(X,Y,VSS); title(iteration); pause(dt);
end
U = US; V = VS;
t = t+dt;
end