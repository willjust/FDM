clf; 
N = 10000;
xx = linspace(0,pi);
y = fpp(xx);
plot(xx,y);
%not sure how to solve a boudary value problem in MATLAB
%sol = bvp4c(@fpp, [0 pi], [0 0]);


function y = fpp(x)
y = 0;
n = 10000;

for k = 1:n
    y = y + 1/k*sin(k*pi/4)^2 * sin(2*k*x);
end
y = 4*y/pi;
end