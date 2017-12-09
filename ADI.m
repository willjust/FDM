function [ T ] = ADI(r, Tinit, time, dt, kappa)
%ADI on initial values given diffusion Kappa
%   Performs ADI on initial values given a diffusion constant kappa
%   Created 10/30/2017 by Justin Williams
x = r(1,:); y = r(2,:); dx = x(2)-x(1); dy = y(2)- y(1);
nx = size(x,2); ny = size(y,2);
alpha_x = kappa*dt/dx^2; alpha_y = kappa*dt/dy^2;

%*************************************************************************%
%                      A matrix formations                                %
%*************************************************************************%
e = ones(nx-2,1);
Ax = spdiags([-alpha_x*e 2*e*(1+alpha_x) -alpha_x*e], -1:1, nx-2, ny-2);
Ay = spdiags([-alpha_y*e 2*e*(1+alpha_y) -alpha_y*e], -1:1, nx-2, ny-2);

%*************************************************************************%
%                         Simulation loop                                 %
%*************************************************************************%
t = 0;
Tnew1 = zeros(nx,ny);
Told = Tinit;
b = zeros(ny-2,1);
while (t < time)
   %**********************************************************************%
   %                        x-direction - round 1                         %
   %**********************************************************************%
   for i = 2:nx-1
	   for j = 2:ny-1
		% Create b
		b(j-1) = alpha_y*(Told(i,j+1) + Told(i,j-1)) + 2*(1-alpha_y)*Told(i,j);
	   end
	   if t == 0 && i == 2
		  fprintf('(2+ alpha_y del_y^2) f_ij^{n+1/2} for first row: \n');
		  disp(b);
	   end
	   % Get new x
	   Tnew1(i,2:nx-1) = thomas(Ax, b);
   end   
     % figure(6); surf(x,y,Tnew1); title('ADI at each time'); pause;
	 
   %**********************************************************************%
   %                        y-direction - round 1                         %
   %**********************************************************************%
   for i=2:ny-1
       % Create b
	   for j = 2:ny-1
		   b(j-1) = alpha_x*Tnew1(j+1,i) + 2*(1-alpha_x)*Tnew1(j,i) + alpha_x*Tnew1(j-1,i); 
	   end
		if t == 0 && i == 2
		    fprintf('(2+ alpha_x del_x^2) f_ij^{n+1/2} for first row: \n');
			disp(b);
		end
	   % Get new y
	   Told(2:nx-1,i) = thomas(Ay,b);
   end
   %figure(6); surf(x,y,Told); title('ADI at each time'); pause;
   t = t+dt; 
   
   %**********************************************************************%
   %                        y-direction - round 2                         %
   %**********************************************************************%
   for i = 2:ny-1
	   % Create b
	   for j = 2:nx-1
		  b(j-1) = alpha_x*Told(j+1,i) + 2*(1-alpha_x)*Told(j,i) + alpha_x*Told(j-1,i); 
	   end 
	   % Get new y
	   Tnew1(2:nx-1,i) = thomas(Ay,b);
   end   
      %figure(6); surf(x,y,Tnew1); title('ADI at each time'); pause;
	  
   %**********************************************************************%
   %                        x-direction - round 2                                   %
   %**********************************************************************%
   for i=2:nx
	   for j = 2:ny-1
		% Create b
		b(j-1) = alpha_y * Tnew1(i,j+1) + 2*(1-alpha_y)*Tnew1(i,j) + alpha_y*Tnew1(i,j-1);
	   end
	   % Get new x
	   Told(i,2:ny-1) = thomas(Ax, b);
   end
  % figure(6); surf(x,y,Told); title('ADI at each time'); pause;
   t = t+dt; 
end
T = Told;
end

