function [alpha, beta, k, cnt] = bisection_method_with_derivative(A,B, l, f)
    % alpha is ana array that stores a1,a2,...,ai
    % beta the same for b1,b2,...,bi
    % k is the total number of iterations required to finish the algo
    % cnt is the amount of times fi(x) was called to be calculated
    
    % Calculate the minimum n that satisfies the tolerance inequation
    n = ceil(log((B - A) / l) / log(2));
    % Derivative of f(x)
    df = diff(f);

    alpha = [];
    beta = [];
    k = 1;
    cnt = 0;

    alpha(end+1) = A;
    beta(end+1) = B;
    for k=1:(n) % k = k + 1 is automatic with `for` loop
        xk = (alpha(end)+beta(end))/2;
        dfxk = df(xk);
        cnt = cnt + 1;
        if dfxk == 0
            break;
        elseif dfxk > 0
            % bhma 2
            alpha(end+1) = alpha(end);
            beta(end+1) = xk;
        else % dfxk < 0
            % bhma 3
            alpha(end+1) = xk;
            beta(end+1) = beta(end);
        end
    end
end

function [] = test(xk, yk, f, grad_f, hessian_f)
    % Tolerance and maximum iterations
    tol = 1e-6;
    max_iter = 500;
    syms g;
    
    % Iterative process
    iter = 0;
    trajectory = [xk, yk]; % To store the trajectory for plotting
    
    % GAMMA WAS CHOSEN THROUGH (B)
    while iter < max_iter
        % Compute gradient and Hessian
        grad = grad_f(xk, yk)';
        hess = hessian_f(xk, yk);
    
        % Check for convergence
        if norm(grad) < tol
            fprintf('Converged to minimum at iteration %d.\n', iter);
            break;
        end
        eign = eig(hess);
        %hess_regularized = hess + 1e-6 * eye(size(hess));
        fprintf("Eigen values of hess are: [%.4f, %.4f]\n", eign(1), eign(2));
        %eig(hess_regularized)
    
        % Check if Hessian is positive definite (to ensure a descent direction)
        %if all(eig(hess_regularized) > 0)
        if all(eig(hess) > 0)
            % Compute Newton direction
            dk = -hess \ grad; % Solves the linear system H * dk = -grad
        else
            % ESSENCE OF THEMA 3
            warning('Hessian is not positive definite. Terminating.');
            break;
        end
        % calculating φ(γ), one variabled function to minimize
        phi(g) = exp(-(xk-g*grad(1))^2 -(yk-g*grad(2))^2)*(xk-g*grad(1))^5;
    
        % calling the minimization technique from Task 1
        % bisection method WITH the use of derivatives
        [a,b,cnt,~] = bisection_method_with_derivative(-10,10, 0.001, phi);
        %fprintf("#%d: cnt=%d\n",k,cnt);
        % finally calculating the gamma that minimizes f(...)
        gamma = (a(end)+b(end))/2;
    
        % Update variables
        xk = xk + gamma * dk(1);
        yk = yk + gamma * dk(2);
        
        % Store the trajectory
        trajectory = [trajectory; xk, yk];
        
        iter = iter + 1;
    end
    
    if iter == max_iter
        fprintf('Reached maximum iterations without convergence.\n');
    end
    
    % Plotting the results
    [X, Y] = meshgrid(-2:0.1:2, -2:0.1:2);
    Z = f(X, Y);
    
    fprintf("Minimum is: [%.6f, %.6f] = %.6f\n", trajectory(end,1), trajectory(end,2), f(trajectory(end,1), trajectory(end,2)));
    figure;
    surf(X, Y, Z);
    hold on;
    plot3(trajectory(:, 1), trajectory(:, 2), f(trajectory(:, 1), trajectory(:, 2)), '-ro', 'LineWidth', 2, 'MarkerSize', 5);
    xlabel('x');
    ylabel('y');
    zlabel('f(x, y)');
    title(sprintf("[%d, %d]: Newton Method Trajectory on Surface", xk(1), yk(1)));
end

clearvars; % clear workspace
clc; % clear command window
close  all;

% Function definition
f = @(x, y) x.^5 .* exp(-x.^2 - y.^2);

% Gradient of f (first derivatives)
grad_f = @(x, y) [5 * x.^4 .* exp(-x.^2 - y.^2) - 2 * x.^6 .* exp(-x.^2 - y.^2), ...
                  -2 * x.^5 .* y .* exp(-x.^2 - y.^2)];

% Hessian matrix (second derivatives)
hessian_f = @(x, y) [20 * x.^3 .* exp(-x.^2 - y.^2) - 22 * x.^5 .* exp(-x.^2 - y.^2) + 4 * x.^7 .* exp(-x.^2 - y.^2), ...
                     4 * x.^6 .* y .* exp(-x.^2 - y.^2) - 10 * x.^4 .* y .* exp(-x.^2 - y.^2); ...
                     4 * x.^6 .* y .* exp(-x.^2 - y.^2) - 10 * x.^4 .* y .* exp(-x.^2 - y.^2), ...
                     4 * x.^5 .* y.^2 .* exp(-x.^2 - y.^2) - 2 * x.^5 .* exp(-x.^2 - y.^2)];


% SOLUTION FOR EVERY INITIAL POINT
test(0,0, f, grad_f, hessian_f);
test(1,-1, f, grad_f, hessian_f);
test(-1,1, f, grad_f, hessian_f);
