% Jonathan Plavnik, CMOR220, Fall 2025, Project 5
% project5_bridge_bonus.m
% Write a function that will solve linear systems based on the LU decomposition
% of a matrix.
% Last modified: October 22, 2025


function [x] = solve_linear_system_lu(A, b)
% list of inputs: A square matrix, b column vector of length 1
% list of outputs: x n by 1 vector solution

[Lower, Upper, Permutation] = lu(A); % LU decomp from PA = LU formula
len = length(b); % size of the col vector 

% Solve Ly = Pb as per the pdf
Pb = Permutation * b; % find Pb
y = zeros(len, 1); % Initialize y w/ length of b
for i = 1:len % We know that L is lower triangular, solve starting at the top 
    sum = 0;

    for j = 1:i-1
        sum = sum + Lower(i, j) * y(j);
    end

    y(i) = (Pb(i) - sum) / (Lower(i, i));
end

% Solve Ux = y
x = zeros(len, 1); % Initializing x with length of b
for i = len:-1:1 % U is upper traingular, solve startin at the bottom 
    sum = 0;
    for j = i+1:len
        sum = sum + Upper(i, j) * x(j);
    end
    x(i) = (y(i) - sum) / (Upper(i, i));
end

end