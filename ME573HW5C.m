%% Problem 1
% ADI method
clc; clf; clear; % I like putting the charts in a specific area, use close all if you  dont;

%% Problem Definition
xMin = 0; xMax = 1; yMin = 0; yMax = 1;
kappa = 0.1; 
dx = 0.05; dy = dx; dt = (dx*dy /kappa)/4; 
t0 = 0; tFinal = 80*dt; 
x = xMin:dx:xMax; y = yMin:dy:yMax; r = [x; y];

%% Solving using different methods
Tinit = tempInit2d(r);
Tftcs = FTCS(Tinit, dx, dy, dt, t0, tFinal, kappa);
Tadi = ADI(r, Tinit, tFinal, dt, kappa);
Tact = tempTest2D(r, kappa, tFinal+dt, 50);
Tact0 = tempTest2D(r, kappa, 0, 50);

%% Visualization/Error
errorSoln = abs(Tact0 - Tinit');
errorADI = abs(Tadi - Tact');
errorFTCS = abs(Tact' - Tftcs);

fprintf('Max Error of Actual = %3.6f\n', norm(reshape(errorSoln, [], 1), inf));
fprintf('Max Error of ADI = %3.6f\n', norm(reshape(errorADI, [], 1), inf));
fprintf('Max Error of FTCS = %3.6f\n', norm(reshape(errorFTCS, [], 1), inf));

figure(1); surf(x,y,Tact'); title('Exact');
figure(2); surf(x,y,errorADI); title('Error');
figure(3); surf(x,y,errorSoln); title('Error of Truncated Solution at T = 0');
hold off