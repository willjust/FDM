%% ME 573 HW 6; Poisson Equation
% By: Justin Williams; November 5, 2017
clc; clf; clear; 

%% Testing
N_test = 2;
A_test= [2 1; 5 7]; b_test = [11; 13];
x_test = jacobi(A_test,b_test);
Error_test = norm(x_test - A_test\b_test);
if Error_test > 1
	fprintf('WARNING: Test case failed for Jacobi.\n');
end


%% Actual Assignment
clc; clf; clear;
xmin = 0; xmax = 1; ymin = 0; ymax = 1; 
dx = 0.01;