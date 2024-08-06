%{
Initialization Funct

Objective: Randomized matrices for A,B,x
     same randomized matrices for all optimization methods
     for compilation comparison purposes
%}

[a,b,x]=init_randn();

function [A,B,x_1] = init_randn()
    dim = input('Enter matrix dimension size dim: ');
    A = randn(dim);
    B = randn(dim,1);
    x_1 = randn(dim,1);

    filename = strcat('randn_matrices_', dim, '.csv');
    
    % csv
    csvwrite(filename, A);
    csvwrite(strcat('B_', filename), B);
    csvwrite(strcat('x1_', filename), x_1);
    
    fprintf('Files created\n');
end
