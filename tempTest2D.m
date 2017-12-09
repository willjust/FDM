function [ Temp ] = tempTest2D(r, kappa, time, trunc)
%TEMPTEST2D approximates the value of x(1-x^5)y(1-y) at time t
%   Created 11/30/2017 by Justin Williams
x = r(1,:); y = r(2,:);
dx = x(2)-x(1); dy = y(2)- y(1);

Tact = zeros(size(x,2),size(y,2));

for n=1:trunc
    for m=1:trunc
        const1 = 120*(-n^4*pi^4*(-1)^n + 12*pi^2*n^2*(-1)^n + 24 + 24*(-1)^(1+n));
        const = -const1*(-2+2*(-1)^m) / (n^7*pi^10*m^3);
		for ii = 1:size(x,2)
            yv = dx*ii-dx;
            for jj=1:size(y,2)
                xv = dy*jj-dy;
                Tact(ii,jj) = Tact(ii,jj) + const*sin(n*pi*xv)*sin(m*pi*yv)*exp(-(n^2+m^2)*pi^2*kappa*(time));
            end
		end
    end
end
Temp = Tact;
end

