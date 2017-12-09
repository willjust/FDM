function [ txy ] = tempInit2d(r)
%TEMPINIT2D Creates a 2D array of temperatures 
%   Creates a 2D temperature profile given by x(1-x^5)y(1-y) given x,y
xx = r(1,:); yy = r(2,:);

nx = size(xx,2); ny = size(yy,2);
dx = xx(2)-xx(1); dy = yy(2)-yy(1);

txy = zeros(nx, ny);

for i=1:nx
	x = i*dx -dx;
	for j=1:ny
		y = j*dy - dy;
		txy(i,j) = x*(1-x^5)*y*(1-y);   
	end
end
end

