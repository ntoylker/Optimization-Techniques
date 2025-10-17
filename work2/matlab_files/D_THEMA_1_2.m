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

clearvars; % clear workspace
clc; % clear command window
close  all;

% Define the function and its gradient
f = @(x, y) x.^5 .* exp(-x.^2 - y.^2);
grad_f = @(x, y) [5*x.^4.*exp(-x.^2-y.^2) - 2*x.^6.*exp(-x.^2-y.^2); ...
                  -2*y.*x.^5.*exp(-x.^2-y.^2)];

% Parameters for the algorithm
initial_points = [0, 0; -1, 1; 1, -1]; % Starting points
tol = 1e-6;                            % Tolerance for convergence
max_iter = 500;                        % Maximum iterations

% PLOTTING f(x,y)
[X, Y] = meshgrid(-2:0.1:2, -2:0.1:2);
Z = f(X, Y);
figure;
surf(X, Y, Z);

% Loop over each starting point
for i = 1:size(initial_points, 1)
    xk = initial_points(i, 1);
    yk = initial_points(i, 2);
    fprintf('\nStarting point: (x0, y0) = (%.1f, %.1f)\n', xk, yk);
    
    % Case (a): Constant step size
    gamma = 0.5; % ADJUST
    trajectory_const = [xk, yk, f(xk, yk)]; % Store trajectory in 3D

    for k = 1:max_iter
        grad = grad_f(xk, yk);
        xk = xk - gamma * grad(1);
        yk = yk - gamma * grad(2);
        trajectory_const = [trajectory_const; xk, yk, f(xk, yk)];
        
        if norm(grad) < tol
            fprintf('Converged with constant step size in %d iterations.\n', k);
            break;
        elseif k==max_iter
            fprintf('Did NOT Converged with constant step size in %d iterations.\n', k);
            break;
        end
    end
    
    % Case (b): minimizing f(xk+γk*dk)  
    % Reset xk and yk for minimization method
    xk = initial_points(i, 1);
    yk = initial_points(i, 2);  
    trajectory_gamma = [xk, yk, f(xk, yk)]; % Store trajectory in 3D
    syms g;

    for k = 1:max_iter
        grad = grad_f(xk, yk);
        % calculating φ(γ), one variabled function to minimize
        phi(g) = exp(-(xk-g*grad(1))^2 -(yk-g*grad(2))^2)*(xk-g*grad(1))^5;

        % calling the minimization technique from Task 1
        % bisection method WITH the use of derivatives
        [a,b,cnt,~] = bisection_method_with_derivative(0,10, 0.001, phi);
        %fprintf("#%d: cnt=%d\n",k,cnt);
        % finally calculating the gamma that minimizes f(...)
        gamma = (a(end)+b(end))/2;

        % xk+1 = ... for next iteration
        xk = xk - gamma * grad(1);
        yk = yk - gamma * grad(2);
        trajectory_gamma = [trajectory_gamma; xk, yk, f(xk, yk)];
        
        if norm(grad) < tol
            fprintf('Converged with step size that minimizes f(xk+γκ*dk) in %d iterations.\n', k);
            break;
        end
    end

    % Case (c): Armijo rule for adaptive step size
    % Reset xk and yk for Armijo method
    xk = initial_points(i, 1);
    yk = initial_points(i, 2);

    trajectory_armijo = [xk, yk, f(xk, yk)]; % Store trajectory for Armijo in 3D
    beta = 0.4;     % Reduction factor
    sigma = 0.01;    % Armijo condition parameter

    for k = 1:max_iter
        grad = grad_f(xk, yk);
        direction = -grad;
        
        % Armijo condition: Find a suitable gamma
        gamma = 1; % Start with an initial guess for gamma

        while f(xk + gamma * direction(1), yk + gamma * direction(2)) > ...
              f(xk, yk) + sigma * gamma * (grad' * direction)
            gamma = beta * gamma; % Reduce gamma if condition is not met
        end
        
        % Update variables
        xk = xk + gamma * direction(1);
        yk = yk + gamma * direction(2);
        trajectory_armijo = [trajectory_armijo; xk, yk, f(xk, yk)];
        
        if norm(grad) < tol
            fprintf('Converged with Armijo step size in %d iterations.\n', k);
            break;
        end
    end
    
    % PLOTTING
    figure;
    % one figure, three plots and a little bit of sauce
    tiledlayout(3,1);
    nexttile;
    k_iter = length(trajectory_const(:,3));
    k = 1:k_iter;
    plot(k, trajectory_const(:,3), '-o');
    title('CONSTANT g', 'Interpreter', 'latex');
    grid on;
    yline(-0.811174, '--', sprintf('-0.811174'), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    xlabel("Iterations");
    ylabel("f(xk)");
    
    nexttile;
    k_iter = length(trajectory_gamma(:,3));
    k = 1:k_iter;
    plot(k, trajectory_gamma(:,3), '-o');
    title("OPTIMAL g", 'Interpreter', 'latex');
    grid on;
    yline(-0.811174, '--', sprintf('-0.811174'), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    xlabel("Iterations");
    ylabel("f(xk)");
    
    nexttile;
    k_iter = length(trajectory_armijo(:,3));
    k = 1:k_iter;
    plot(k, trajectory_armijo(:,3), '-o');
    title("ARMIJO g", 'Interpreter', 'latex');
    grid on;
    yline(-0.811174, '--', sprintf('-0.811174'), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    xlabel("Iterations");
    ylabel("f(xk)");

    % Create a 3D surface plot of the function
    figure;
    x_vals = linspace(-2, 2, 100);
    y_vals = linspace(-2, 2, 100);
    [X, Y] = meshgrid(x_vals, y_vals);
    Z = f(X, Y);
    surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    
    % Plot the trajectory for the constant step size
    plot3(trajectory_const(:,1), trajectory_const(:,2), trajectory_const(:,3), 'r-o', ...
          'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Constant Step');
 
    % Plot the trajectory for the step size that minimizes f(...)
    plot3(trajectory_gamma(:,1), trajectory_gamma(:,2), trajectory_gamma(:,3), 'y-o', ...
          'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Minimization γ Step');

    % Plot the trajectory for the Armijo step size
    plot3(trajectory_armijo(:,1), trajectory_armijo(:,2), trajectory_armijo(:,3), 'b-o', ...
          'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Armijo Step');

    % Add labels and legend
    title(sprintf('3D Trajectories from Starting Point (%.1f, %.1f)', initial_points(i, 1), initial_points(i, 2)));
    xlabel('x');
    ylabel('y');
    zlabel('f(x, y)');
    legend;
    grid on;
    hold off;
end
