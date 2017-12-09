function y = fcn(x, n, t, L, k) 
PI = 3.14159;
y = 0;
for i=1:n
    cn = 4*L / (PI^2 * i^2)*sin(PI*i/2);
    y = y + cn*sin(i*PI*x/L)*exp(-(i*PI*sqrt(k)/L)^2 * t); 
end