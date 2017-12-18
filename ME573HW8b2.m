clear; clf;

%% Setup domain and get exact solution
xmin = -1; xmax = 1; ymin = 0; ymax = 2;
dx = 0.1; dy = dx; t0 = 0; dt = 0.025; tf = 5;
gx = xmin:dx:xmax; gy = ymin:dy:ymax;
[X, Y] = meshgrid(gx, gy);
[U, V] = FletchersExact(X,Y,-1);
figure(1); surf(X,Y,U); xlabel('x'); ylabel('y'); title('Exact U')
figure(2); surf(X,Y,V); xlabel('x'); ylabel('y'); title('Exact V')

kappa = 1/5;
%% FTCS Solution
[Uf, Vf] = FTCS(zeros(size(X)), X, Y, [dx, dy], dt, tf, kappa, @FletchersExact, 'nonlinear');
figure(3); surf(X,Y,abs(U-Uf)); xlabel('x'); ylabel('y'); title('FTCS Error-U')
figure(4); surf(X,Y,abs(V-Vf)); xlabel('x'); ylabel('y'); title('FTCS Error-V')
printc('FTCS Max Error');
fue = norm(reshape(U-Uf, [], 1), inf);
fve = norm(reshape(V-Vf, [], 1),inf);
fprintf('U: %3.3f\tV: %3.3f\n', fue, fve);

%% Crank-Nicolson Solution
Uc = U; Vc = V; nx = size(X,1); ny = size(X,2);
Uc(2:nx-1, 2:ny-1) = zeros(nx-2, ny-2);
Vc(2:nx-1, 2:ny-1) = zeros(nx-2, ny-2);
[Uc, Vc] = crankNicolson(Uc,Vc,dx,dy,dt,tf,kappa);
%[Uc, Vc] = crankNicolsen(zeros(size(X)), X, Y, [dx, dy], dt, tf, kappa, @FletchersExact);
figure(5); surf(X,Y,abs(U-Uc)); xlabel('x'); ylabel('y'); title('CN Error-U')
figure(6); surf(X,Y,abs(V-Vc)); xlabel('x'); ylabel('y'); title('CN Error-V')

printc('Crank-Nicolson Results');
cue = norm(reshape(U-Uc, [], 1), inf);
cve = norm(reshape(V-Vc, [], 1), inf);
fprintf('U: %3.3f\tV: %3.3f\n', cue, cve);
