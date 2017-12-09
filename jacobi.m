function [ x ] = jacobi( A, b, varargin)
%Jacobi performs an iterative solver on a matrix
%   Performs Jocobi iteration method on a matrix using
%	x^{(k+1)} = Tx^{(k)} + C
%	Created: Justin Williams, November 5, 2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%						Setup Parameters								  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(A, 1);

% Setting up Equation x^{(k+1)} = Tx^{(k)} + C, where T and C are const.
T = zeros(n,1);
for i=1:n
	T(i,i) = 1/A(i,i); % D^{-1}
	A(i,i) = 0; % Make A = L+U
end
C = T*b;  % C = D^{-1}*b (T is still D  here)
T = -T*A; % T = D^{-1} (L+U)

iterations = 50;
if ~isempty(varargin)
	iterations = varargin{1};
end
x = zeros(n,1);
error = norm(b-x); % inital L_2 norm error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%						Iteratively Solve								  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:iterations
	xp = x;
	x = T*x + C; 
	error = norm(x-xp);
end
if error > 10
	fprintf('WARNNING: Jacobi not converged, residual: %d\n', error);
end
end

