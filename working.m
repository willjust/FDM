L = 1;
k = 0.001;
N = 200;
t_med = 0.1*L^2 / (pi^2*k);
t_end=0.5*L^2 / (pi^2*k);
xx = linspace(0,L);
tt = [0, t_med, t_end];

figure(1)
title('function value at various t');
xlabel('x');
ylabel('fxn(x,t)');
hold on
yy1 = fcn(xx, N, tt(1), L, k);
yy2 = fcn(xx, N, tt(2), L, k);
yy3 = fcn(xx, N, tt(3), L, k);
plot(xx,yy1,xx, yy2,xx,yy3);
hold off

figure(2)
hold on
title('Error vs N');
xlabel('N');
ylabel('Error');
yy_1_exact = fcn(xx, 10000, tt(1), L, k);
yy_2_exact = fcn(xx, 10000, tt(2), L, k);
yy_3_exact = fcn(xx, 10000, tt(3), L, k);

err1 = zeros(50);
err2 = zeros(50);
err3 = zeros(50);

for n = 1:50
yy11 = fcn(xx, n, tt(1), L, k);
yy21 = fcn(xx, n, tt(2), L, k);
yy31 = fcn(xx, n, tt(3), L, k);

err1(n) = norm(yy_1_exact - yy11, inf);
err2(n) = norm(yy_2_exact - yy21, inf);
err3(n) = norm(yy_3_exact - yy31, inf);
end
n = 1:50;
plot(n, err1, n, err2, n, err3);