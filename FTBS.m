function y = FTBS(Finit, dx, dt, tf, kappa, U)
	dim = size(Finit, 1);
	y = Finit;
	if(dim==1)
		y = FTBS1D(Finit, dx, dt, tf, kappa, U);
	elseif (dim==2)
		y = FTBS2D(Finit,dx,dt,tf,kappa,U); % not actually implemented yet
	elseif (dim == 3)
		fprintf('ERROR: 3D FTBS is not implemented.\n');
	end
end

function y = FTBS1D(Finit, dx, dt, tf, kappa, U) 
	% Check upwindiness of method
	negate = 0;
	if (U < 0)
		U = -U;
		negate = 1;
	end
	n = size(Finit,2); F1 = Finit; F2=Finit; t = 0;
	alpha = kappa*dt/dx^2; CFL = U*dt/dx;
	if(CFL^2 > CFL + 2*alpha || CFL + 2*alpha > 1)
		fprintf('WARNING: FTBS solution unstable! Use different method.\n')
		fprintf('FTBS -> CFL = %d \t alpha = %d\n', CFL, alpha);
	end
	
	while(t<tf)
		for i = 2:n-1
			F2(i) = F1(i-1) * (CFL+alpha) + F1(i)*(1-CFL - 2*alpha) ...
				+ F1(i+1)*alpha;
		end
		F1 = F2; t = t+dt;
	end
	y = F1;
	
	if(negate)
		for i = 1:n
			y(i) = F1(n-i+1);
		end
	end
end

% TODO
function y = FTBS2D(Finit, dx, dt, tf, kappa, U)
% FTBS Performs 2D FTBS scheme on Finit
%	Created 10/26/2017 by Justin Williams
	if(U ~= 0)
		fprintf('ERROR: Convection not implemented in 2D FTBS\n');
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
        fprintf('WARNING: FTBS in unstable regime\n');
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
    fprintf('2D FTBS Finished.\n');
    y = T1;
end
