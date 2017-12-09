%% Thomas' algorithm for solving a tridiagonal matrix
% Created by Justin Williams
function [y] = thomasAlg(A, b) 
    N = size(A,1);
    l = zeros(N);
    for k=1:N-1
        i = k+1;
        l(i,k) = A(i,k) / A(k,k);
        for j = k:k+1
            A(i,j) = A(i,j) - l(i,k)*A(k,j);
        end
        b(i) = b(i) - l(i,k)*b(k);
    end
    y = zeros(N,1);
    for k=N:-1:1
        y(k) = b(k);
        for j=k+1:min(N,k+1)
            y(k) = y(k) - A(k,j)*y(j);
        end
        y(k)=y(k)/A(k,k);
    end
end