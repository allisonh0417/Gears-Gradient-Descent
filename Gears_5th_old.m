%{
Gears method (1st->5th Order) with Adaptive Step Size
Backward Differentiation 
equations provided by https://www.chemecomp.com/Gear.pdf 

Directions: Run 'init_func.m' first to generate random NxN matrices
%}

A = a;
B = b;
y_curr = x; % y(n+1)

k = 0; % iteration counter
z = 0; % switch; if z =0, keep increasing step size; if z=1, decrease
epsilon = 1e-6; %value close to 0 for convergence
num = 1; %(num)th Order 
last_valid_err = 0; %last valid error in function in case need to go back to 1st Order

h_1 =1/(3*norm(A)*norm(A)); % step size h_1

g = f(y_curr, A, B);  %gradient function
err = norm(A * y_curr - B); %error function

%initialize plotting vectors
error_values = [];
norm_values = [];

tic %start timer
while err > epsilon  
    if mod(k, 10000) == 1 
            fprintf('Iteration: %d, Error: %.10f, num: %d, z:%d, Step size: %.10f\n', k, err, num, z,h_1);
    end
    k = k + 1;

    if num == 1 %1st Order
        y_prev_0 = y_curr; % y(n+1) from prev. iteration -> y(n)
        y_curr_test = y_prev_0 + h_1 * f(y_prev_0, A, B); % new (potential) y(n+1) with adjusted h_1
        err_1 = norm(A * y_curr_test - B); % new (potential) error value
        if err_1 >= err %verify if new error is smaller than previous iteration's error
            h_1 = h_1*0.9; %if too big, decrease step size
        else
            num = 2;
        end
        y_curr = y_curr_test;
        last_valid_err = err_1;
    elseif num == 2 %2nd Order
        y_prev_1 = y_curr; 
        y_curr = (4/3)*y_prev_1 - (1/3)*y_prev_0 + (2/3)*h_1*f(y_prev_1, A, B);
        err_2 = norm(A * y_curr - B);
        if err_2 > last_valid_err 
            num = 1; 
            [h_1, z] = new_step(z, err, h_1, A, B, y_curr);
        else
            last_valid_err = err_2;
            num = 3;
        end
    elseif num == 3 %3rd Order
        y_prev_2 = y_curr; 
        y_curr = (18/11)*y_prev_2 - (9/11)*y_prev_1 + (2/11)*y_prev_0 + (6/11)*h_1*f(y_prev_2, A, B);
        err_3 = norm(A * y_curr - B);
        if err_3 > last_valid_err 
            num = 1; 
            [h_1, z] = new_step(z, err, h_1, A, B, y_curr);
        else
            last_valid_err = err_3;
            num = 4;
        end
    elseif num == 4 %4th Order
        y_prev_3 = y_curr;
        y_curr = (48/25)*y_prev_3 - (36/25)*y_prev_2 + (16/25)*y_prev_1 - (3/25)*y_prev_0 + (12/25)*h_1*f(y_prev_3, A, B);
        err_4 = norm(A * y_curr - B);
        if err_4 > last_valid_err 
            num = 1; 
            [h_1, z] = new_step(z, err, h_1, A, B, y_curr);
        else
            last_valid_err = err_4;
            num = 5;
        end
    else %5th Order
        y_prev_4 = y_curr; %is this okay?
        y_curr = (300/137)*y_prev_4 - (300/137)*y_prev_3 + (200/137)*y_prev_2 - (75/137)*y_prev_1 + (12/137)*y_prev_0 + (60/137)*h_1*f(y_prev_4, A, B);
        err_5 = norm(A * y_curr - B);
        if err_5 > last_valid_err 
            [h_1, z] = new_step(z, err, h_1, A, B, y_curr);
            num = 1; 
        else
            last_valid_err = err_5;
        end
    end

    % Plot values
    error_values = [error_values; err];
    norm_values = [norm_values; norm(y_curr)];

    err = norm(A * y_curr - B);
end

%y_curr,g
k
norm_x = norm(y_curr)
err = norm(A*y_curr-B)
toc

%
% Plot log[error] values
figure;
semilogy(error_values, 'LineWidth', 1.5); % Use semilogy for logarithmic plot
xlabel('Iteration');
ylabel('Log(Error)');
title('Error Convergence (Log Scale)');
grid on;

% Plot norm[y_curr] values, aka norm[X]
figure;
plot(norm_values, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Norm(y\_curr)');
title('Norm of y\_curr Convergence');
grid on;

%}

function [f_val] = f(y, A, B) %gradient function
    f_val = -2 * A' * (A * y - B);
end

function [h_new, z] = new_step(z, err, h_1, A, B, y_curr) %adaptive step size function
    if z == 0
        h_test = h_1 * 1.1;
        y_new_test = y_curr + h_test * f(y_curr, A, B);
        err_new = norm(A * y_new_test - B);
        if err_new > err 
            z = 1;
        else
            h_1 = h_test; 
        end
    else
        h_test = h_1 * 0.95;
        y_new_test = y_curr + h_test * f(y_curr, A, B);
        err_new = norm(A * y_new_test - B);
        if err_new > err
            z = 0;
        else
            h_1 = h_test; %only change step size if new estimated error is smaller than current error
        end
    end
    h_new = h_1;
end
