%{
Gears method (1st->5th Order) with Adaptive Step Size
Backward Differentiation 
equations provided by https://www.chemecomp.com/Gear.pdf 

Directions: Run 'init_func.m' first to generate random NxN matrices
%}
%This one will always have decreasing error

A = a;
B = b;
y_curr = x; 

k = 1; % iteration counter
z = 0; % switch; if z =0, keep increasing step size; if z=1, decrease
epsilon = 1e-6; %value close to 0 for convergence
num = 1; %(num)th Order 
g=0; %counter in step func
clk=0; %global counter, clk=k+g

h_1 =1/(3*norm(A)*norm(A)); % step size h_1

grad = f(y_curr, A, B);  %gradient function
err = norm(A * y_curr - B); %error function
last_valid_err = err; %last valid error in function in case need to go back to 1st Order


%initialize plotting vectors
error_values = [err];
norm_values = [norm(y_curr)];
step_size_values = [h_1];


tic %start timer
while err > epsilon 
    t=t+1;
    if num == 1 %1st Order
        y_prev_0 = y_curr; % y(n+1) from prev. iteration -> y(n)
        y_curr_test = y_prev_0 + h_1 * f(y_prev_0, A, B); % new (potential) y(n+1) with adjusted h_1
        err_1 = norm(A * y_curr_test - B); % new (potential) error value
        if err_1 >= last_valid_err %verify if new error is smaller than previous iteration's error
                [h_1, z,g] = new_step(z, err_1, last_valid_err,h_1, A, B, y_curr);
                num=2;
                y_curr = y_curr_test;
        else
            num = 2;
            y_curr = y_curr_test;
            last_valid_err = err_1;
            k=k+1;
        end
        
    elseif num == 2 %2nd Order
        y_prev_1 = y_curr; 
        y_curr = (4/3)*y_prev_1 - (1/3)*y_prev_0 + (2/3)*h_1*f(y_prev_1, A, B);
        err_2 = norm(A * y_curr - B);
        if err_2 > last_valid_err 
            num = 1; 
            [h_1, z,g] = new_step(z, err_2, last_valid_err,h_1, A, B, y_curr);
        else
            last_valid_err = err_2;
            num = 3;
            k=k+1;
        end
    elseif num == 3 %3rd Order
        y_prev_2 = y_curr; 
        y_curr = (18/11)*y_prev_2 - (9/11)*y_prev_1 + (2/11)*y_prev_0 + (6/11)*h_1*f(y_prev_2, A, B);
        err_3 = norm(A * y_curr - B);
        if err_3 > last_valid_err 
            num = 1; 
            [h_1, z,g] = new_step(z, err_3,last_valid_err, h_1, A, B, y_curr);
        else
            last_valid_err = err_3;
            num = 4;
            k=k+1;
        end
    elseif num == 4 %4th Order
        y_prev_3 = y_curr;
        y_curr = (48/25)*y_prev_3 - (36/25)*y_prev_2 + (16/25)*y_prev_1 - (3/25)*y_prev_0 + (12/25)*h_1*f(y_prev_3, A, B);
        err_4 = norm(A * y_curr - B);
        if err_4 > last_valid_err 
            num = 1; 
            [h_1, z,g] = new_step(z, err_4,last_valid_err, h_1, A, B, y_curr);
        else
            last_valid_err = err_4;
            num = 5;
            k=k+1;
        end
    else %5th Order
        y_prev_4 = y_curr; %is this okay?
        y_curr = (300/137)*y_prev_4 - (300/137)*y_prev_3 + (200/137)*y_prev_2 - (75/137)*y_prev_1 + (12/137)*y_prev_0 + (60/137)*h_1*f(y_prev_4, A, B);
        err_5 = norm(A * y_curr - B);
        if err_5 > last_valid_err 
            [h_1, z,g] = new_step(z, err_5,last_valid_err, h_1, A, B, y_curr);
            y_curr = (300/137)*y_prev_4 - (300/137)*y_prev_3 + (200/137)*y_prev_2 - (75/137)*y_prev_1 + (12/137)*y_prev_0 + (60/137)*h_1*f(y_prev_4, A, B);
            num = 1;
        else
            last_valid_err = err_5;
            k=k+1;
        end
    end
    if mod(k, 1000) == 1  
            fprintf('Iteration: %d, Error: %.10f, num: %d, z:%d, LastV Error: %.10f, ctr: %d,Step size: %.10f\n', k, err, num, z, last_valid_err,clk,h_1);
    end

    err = norm(A * y_curr - B);
    % Plot values
    if num ~=1 &&last_valid_err >= err
        error_values = [error_values; last_valid_err];
        norm_values = [norm_values; norm(y_curr)];
        step_size_values = [step_size_values;h_1];
    end

    clk = g +clk +1;
end

% Results
norm(y_curr)
err
k
toc

gear_error = error_values;
gear_norm = norm_values;

%{
figure; 

% Plot step_size values on the right y-axis first
yyaxis left;
lighter_red = [1, 0.6, 0.6]; % Define a lighter shade of red
plot(step_size_values,  'LineWidth', 1.5, 'Color', lighter_red); % Plot step_size with lighter red
ylabel('Step Size', 'Color', lighter_red); % Label in the same lighter red

% Plot log[error] values on the left y-axis, after the step_size plot
yyaxis right;
semilogy(gear_error, 'LineWidth', 1.5, 'Color', 'b'); % Plot log[error] in blue
xlabel('Iteration');
ylabel('Log(Error)', 'Color', 'b'); % Label in blue
title('Gears NEW: Error and Step Size Convergence');
grid on;

% Add a legend to identify the plots
legend({'Step Size', 'Log(Error)'}, 'Location', 'best');

% Adjust the y-axis properties for readability
ax = gca; % Get current axis
ax.YAxis(1).Color = lighter_red; % Set the left y-axis color to blue
ax.YAxis(2).Color = 'b'; % Set the right y-axis color to the lighter red
%}

%{
% Plot log[error] values
figure;
semilogy(gear_error, 'LineWidth', 1.5); % Use semilogy for logarithmic plot
xlabel('Iteration');
ylabel('Log(Error)');
title('Error Convergence (Log Scale): Gears NEW');
grid on;
%

% Plot norm[y_curr] values, aka norm[X]
figure;
plot(gear_norm, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Norm(y\_curr)');
title('Norm of y\_curr Convergence: Gears NEW');
grid on;

% Plot step_size values
figure;
plot(step_size_values, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('step size');
title('Step Size Value: Gears NEW');
grid on;
%


% Plot norm[y_curr] values, aka norm[X]
figure;
plot(gear_norm, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Norm(y\_curr)');
title('Norm of y\_curr Convergence: Gears NEW');
grid on;
%}

function [f_val] = f(y, A, B) %gradient function
    f_val = -2 * A' * (A * y - B);
end

%
function [h_new, z,t] = new_step(z, local_err, global_err,h_1, A, B, y_curr) %adaptive step size function
    t=1;
    err_new = local_err;
    err_prev = local_err;
    while err_new >= global_err  && t<5
        
        if z == 0 %increase
            h_test = h_1 * 1.1;
            y_new_test = y_curr + h_test * f(y_curr, A, B);
            err_new = norm(A * y_new_test - B);
            if err_new < local_err
                local_err = err_new;
                %t=t+1;
                h_1 = h_test;
                err_prev=err_new;
            elseif err_new < global_err
                h_1 = h_test; 
                break
            elseif err_prev < err_new
                t=t+1;
                z=1;
                err_prev=err_new;
            else
                t=t+1;
                err_prev=err_new;
            end
        else %decrease
            h_test = h_1 * 0.9;
            y_new_test = y_curr + h_test * f(y_curr, A, B);
            err_new = norm(A * y_new_test - B);
            if err_new < local_err
                local_err = err_new;
                h_1 = h_test;
                err_prev=err_new;
            elseif err_new<global_err
                h_1 = h_test; 
                break
            elseif err_prev < err_new
                t=t+1;
                z=0;
                err_prev=err_new;
            else
                t=t+1;
                err_prev=err_new;
            end
        end
    end
    h_new = h_1;
end


