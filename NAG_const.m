%Nesterov's Adaptive Gradient 
%constant step size
%Objective: try to find most optimal step size
% NOTE: Initialize A, B, and x before running the script

A = a;
B = b;
y_curr = x;

k = 1;
epsilon = 1e-6; 

h_1=0.001; %constant step size, insert (if too big will diverge)
%Tip: you can use NAG_Adaptive to find range for optimal step size h
err = norm(A * y_curr - B);

step_values = [h_1];
error_values = [err];
norm_values = [norm(y_curr)];
tic
while err > epsilon 
    if mod(k, 10000) == 1
        fprintf('Iteration: %d, Error: %.10f, Step size: %.10f\n', k, err, h_1);
    end
    
    k =k+1;
    %NAG
    if k==1
        y_prev = y_curr;
    end
    y_curr = y_prev - h_1 *f(y_curr,A,B);
    y_new = y_curr + ((k-1)/(k+2))*(y_curr-y_prev);
    
    y_prev = y_curr;
    y_curr = y_new;

    err = norm(A * y_curr - B);

    step_values = [step_values; h_1];
    error_values = [error_values; err];
    norm_values = [norm_values; norm(y_curr)];
end

% Final results
y_curr, g;
k
err = norm(A * y_curr - B)
toc

NAG_step_val = step_values;
NAG_error = error_values;
NAG_norm = norm_values;

function [f_val] = f(y, A, B)
    f_val = 2 * A' * (A * y - B);
end
