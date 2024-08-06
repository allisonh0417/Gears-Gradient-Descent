%Nesterov's Adaptive Gradient 
%step size w adaptive adjustment

% NOTE: Initialize A, B, and x before running the script

A = a;
B = b;
y_curr = x;
z = 0;

k = 0;
epsilon = 1e-6; 
tic 

%h_1 = 0.001; % Initial step size
h_1=1/(3*norm(A)*norm(A));
err = norm(A * y_curr - B);
g = f(y_curr, A, B);

error_values = [];
norm_values = [];

while err > epsilon
    if mod(k, 10000) == 0
        fprintf('Iteration: %d, Error: %.10f, Step size: %.10f\n', k, err, h_1);
    end

    k = k + 1;

    %NAG
    y_prev = y_curr;
    y_curr = y_prev - h_1 *f(y_curr,A,B);
    y_new = y_curr + ((k-1)/(k+2))*(y_curr-y_prev);
    g = f(y_curr, A, B);

    err = norm(A * y_curr - B);

    error_values = [error_values; err];
    norm_values = [norm_values; norm(y_curr)];

    % Adaptive step size adjustment
    if z == 0
        h_test = h_1 * 1.005; % Increase step size
    else
        h_test = h_1 * 0.9; % Decrease step size
    end

    %NAG
    y_prev_test = y_curr;
    y_curr_test = y_prev_test - h_test *f(y_curr,A,B);
    y_new_test = y_curr_test + ((k-1)/(k+2))*(y_curr_test-y_prev_test);

    err_new = norm(A * y_new_test - B);

    if err_new < err
        % Accept the new step and update h_1
        h_1 = h_test;
        err = err_new;
        y_curr = y_new_test;
        z = 0; % Continue increasing step size
    else
        % Revert to old step size
        z = 1; % Start decreasing step size
    end
end

% Final results
y_curr, g, k
err = norm(A * y_curr - B)
toc

% Plotting the error values with logarithm transformation
figure;
semilogy(error_values, 'LineWidth', 1.5); % Use semilogy for logarithmic plot
xlabel('Iteration');
ylabel('Log(Error)');
title('Error Convergence (Log Scale)');
grid on;

% Plotting the norm of y_curr values
figure;
plot(norm_values, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Norm(y\_curr)');
title('Norm of y\_curr Convergence');
grid on;

function [f_val] = f(y, A, B)
    f_val = 2 * A' * (A * y - B);
end