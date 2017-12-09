function y = LinAdvecAnalytical(x, c0, ci, D, R, v, t)
A = @(x,t) 1/2*erfc((R*x-v*t)/(2*sqrt(D*R*t))) + 1/2 * (exp(v*x/D)*erfc((R*x+v*t)/(2*sqrt(D*R*t))));
y = zeros(size(x));
for i = 1:size(x,2)
	y(i) = ci +	(c0-ci)*A(x(i),t);
end
end