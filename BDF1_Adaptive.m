%Gears method Adaptive (adaptive step-size)
%Backward Differentiation Method

A = a;
B = b;
y_curr = x;
k=0;
z=0; % if z =0, keep increasing step size
epsilon = 0.000001;
val=0;

%func1 = (transpose(x)*transpose(A)*A*x) - 2*transpose(B)*A*x;
%grad1 = 2*transpose(A)*A*x - 2*transpose(A)*B;
tic 

%h_1 = 0.001; %step size, currently constant
h_1=1/(3*norm(A)*norm(A));
g=f(y_curr,A,B);  

initial_err = norm(A * y_curr - B);
error_values = []; 
norm_values = []

while norm(A * y_curr - B) > epsilon 
    %if (mod(k,10000)==1) fprintf('%f %.5f \n',k,y_curr); end
    if mod(k, 10000) == 1
        fprintf('Iteration: %d, Error: %.10f, val: %d,Step size: %.10f\n', k, err, val,h_1);
    end
    k = k + 1;
    y_prev = y_curr; % y_curr is old value -> prev.
    y_curr = y_prev + (h_1) * f(y_curr, A, B);

    % adaptive step size h
    err = norm(A * y_curr - B);
   
    if err<= initial_err
         error_values = [error_values; err]; % Store current error value
         norm_values = [norm_values; norm(y_curr)]; % Store the norm of y_curr
    val = norm(y_curr);
    end

    if z == 0
        if err > epsilon
            h_inc = h_1*1.1;%inc step size
            h_test = h_inc; 
            % test out new step size, check if error increased
            y_prev_test = y_curr; % y_curr is old value -> prev.
            y_new_test = y_prev + (h_test) * f(y_prev, A, B);
            err_new = norm(A * y_new_test - B);

            if err_new > err  % revert back to old values if error increases
                z = 1;
            else
                h_1 = h_test; % update step size if error does not increase
            end
        end
            
    else %need to start decreasing step size
        if err > epsilon
            h_inc = h_1*0.9; %decrease step size
            h_test = h_inc;
            % test out new step size, check if error increased
            y_prev_test = y_curr; % y_curr is old value -> prev.
            y_new_test = y_prev + (h_test) * f(y_prev, A, B);
            err_new = norm(A * y_new_test - B);

            if err_new > err  % revert back to old values if error increases
                z = 0;
            else
                h_1 = h_test; % update step size if error does not increase
            end
        end       
    end
    g = f(y_curr, A, B);
end



y_curr,g,k
err = norm(A*y_curr-B)
toc

% Plotting the error values
figure;
semilogy(error_values, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Error');
title('Error Convergence');
grid on;

figure;
subplot(2,1,2);
plot(norm_values, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Norm(y\_curr)');
title('Norm of y\_curr Convergence');
grid on;


function [f_val] = f(y,A,B) %gradient descent funct.
    f_val = -2*A'*(A*y-B);
end