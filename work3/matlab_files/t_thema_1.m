function [trajectory] = steepest_descent(f, grad_f, x0, tol, gamma, max_iter)

    grad = grad_f(x0(1), x0(2));
    x_arr = x0(1);
    y_arr = x0(2);
    trajectory = [x0(1), x0(2), f(x0(1), x0(2))]; % Store trajectory in 3D
    k=0;

    
    % MAIN ALGO IMPLEMENTATION

    while norm(grad)>tol
        grad = grad_f(x_arr, y_arr);
        x_arr = x_arr - gamma * grad(1);
        y_arr = y_arr - gamma * grad(2);
        trajectory = [trajectory; x_arr, y_arr, f(x_arr, y_arr)];
        
        k=k+1;
        if k>max_iter
            fprintf('Did NOT Converge with constant step size g=%.3f in %d iterations.\n', gamma,k-1);
            break;
        end 
    end

    if k<max_iter && trajectory(end, 3) < 0.1
        fprintf('Converged with constant step size g=%.3f in %d iterations.\n', gamma, k-1);
    elseif k<max_iter
        fprintf('Did NOT converge. Went to infinity for constant step sizes g=%.3f in %d iterations.\n', gamma, k-1);
    end



    % PLOTTING
    figure;
    % Create a 2D contour plot
    if x0(1) > 0
        x_vals = linspace(-x0(1), x0(1), 100);
    else
        x_vals = linspace(x0(1), -x0(1), 100);
    end

    if x0(2) > 0
        y_vals = linspace(-x0(2), x0(2), 100);
    else
        y_vals = linspace(x0(2), -x0(2), 100);
    end

    [X, Y] = meshgrid(x_vals, y_vals);
    Z = f(X, Y);
    contour(X, Y, Z, 50, 'LineWidth', 1.2);
    hold on;
    % Plot the trajectory on the contour plot
    plot(trajectory(:, 1), trajectory(:, 2), 'r-o', 'LineWidth', 2, 'MarkerSize', 5, ...
         'DisplayName', 'Trajectory');
    scatter(trajectory(1, 1), trajectory(1, 2), 100, 'g', 'filled', 'DisplayName', 'Start Point');
    scatter(trajectory(end, 1), trajectory(end, 2), 100, 'b', 'filled', 'DisplayName', 'End Point');
    % Add labels and legend
    title(sprintf('Contour Plot with Trajectory (gamma = %.3f)', gamma), 'Interpreter', 'latex');
    xlabel('x', 'Interpreter', 'latex');
    ylabel('y', 'Interpreter', 'latex');
    legend('show', 'Interpreter', 'latex');
    grid on;



    figure;
    k_iter = length(trajectory(:,3));
    k = 1:k_iter;
    plot(k, trajectory(:,3), '-o');
    title(sprintf('gamma = %.3f',gamma), 'Interpreter', 'latex');
    grid on;
    %yline(-0.811174, '--', sprintf('-0.811174'), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    xlabel("Iterations");
    ylabel("f(xk)");
    


    figure;
    % Create a 3D surface plot of the function
    if x0(1) > 0
        x_vals = linspace(-x0(1), x0(1), 100);
    else
        x_vals = linspace(x0(1), -x0(1), 100);
    end

    if x0(2) > 0
        y_vals = linspace(-x0(2), x0(2), 100);
    else
        y_vals = linspace(x0(2), -x0(2), 100);
    end
    
    [X, Y] = meshgrid(x_vals, y_vals);
    Z = f(X, Y);
    surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    hold on;    
    % Plot the trajectory for the step size
    plot3(trajectory(:,1), trajectory(:,2), trajectory(:,3), 'r-o', ...
          'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Trajectory');
    % Add labels and legend
    title(sprintf('3D Trajectories from Starting Point (%.1f, %.1f)', x0(1), x0(2)));
    xlabel('x');
    ylabel('y');
    zlabel('f(x, y)');
    legend;
    grid on;



    figure;
    k_iter = length(trajectory(:,3));
    k = 1:k_iter;
    plot(k, trajectory(:,1));
    hold on;
    plot(k, trajectory(:,2));

    title(sprintf('gamma = %.3f',gamma), 'Interpreter', 'latex');
    grid on;
    xlabel("Iterations k");
    ylabel("x, y Common Axes");
    legend('x', 'y'); % Set legend labels
end




clearvars; % clear workspace
clc; % clear command window
close  all;

% Define the function and its gradient
f = @(x, y) 1/3*x.^2 + 3*y.^2;
grad_f = @(x, y) [2/3*x; 6*y];

% PLOTTING f(x,y)
[X, Y] = meshgrid(-2:0.1:2, -2:0.1:2);
Z = f(X, Y);
figure;
surf(X, Y, Z);

% Parameters for the algorithm
initial_points = [4, 4]; % Starting points
tol = 1e-3;                             % Tolerance for convergence
max_iter = 2000;                        % Maximum iterations

gamma =  [0.1, 0.3, 1/3-tol, 1/3+tol, 3, 5];

steepest_descent(f, grad_f, initial_points(1,:), tol, gamma(1), max_iter);
steepest_descent(f, grad_f, initial_points(1,:), tol, gamma(2), max_iter);
steepest_descent(f, grad_f, initial_points(1,:), tol, gamma(3), max_iter);
steepest_descent(f, grad_f, initial_points(1,:), tol, gamma(4), max_iter);
steepest_descent(f, grad_f, initial_points(1,:), tol, gamma(5), max_iter);
steepest_descent(f, grad_f, initial_points(1,:), tol, gamma(6), max_iter);
