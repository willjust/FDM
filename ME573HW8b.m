%% ME573 HW8b.
%Justin Williams 12/6/2017
% Problem 2
dx = 0.1; dy = 0.1; dt = 0.001;
xmin = -1; xmax = 1; x = xmin:dx:xmax; nx = size(x,1);
ymin = 0; ymax = 2; y = ymin:dy:ymax; ny = size(y,1);
tf = 10;
k = 1/2; nu = 1/5; x0 = 1; a0 = 0.001*k*exp((1+x0)*k); a1 = a0; 
ymax = ymax-ymin; a2 = 0; a3 = 0;

