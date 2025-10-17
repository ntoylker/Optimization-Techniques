clearvars; % clear workspace
clc; % clear command window
close  all;

f = @(x, y) x.^5 .* exp(-x.^2 - y.^2);

grad_f = @(x, y) [5*x.^4.*exp(-x.^2-y.^2) - 2*x.^6.*exp(-x.^2-y.^2); ...
                  -2*y.*x.^5.*exp(-x.^2-y.^2)];
hessian_f = @(x, y) [ 20*x^3*exp(- x^2 - y^2) - 22*x^5*exp(- x^2 - y^2) + ...
                          4*x^7*exp(- x^2 - y^2), 4*x^6*y*exp(- x^2 - y^2) - 10*x^4*y*exp(- x^2 - y^2) ;
                        4*x^6*y*exp(- x^2 - y^2) - 10*x^4*y*exp(- x^2 - y^2),  4*x^5*y^2*exp(- x^2 - y^2) - 2*x^5*exp(- x^2 - y^2) ];


function [alpha, beta, k, cnt] = bisection_method_derivatives(A, B, l, f)
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


function [x_vals, y_vals, f_vals] = lm_method_fixed(  f, grad_f, hessian_f, x0, epsilon, gamma)
    x = x0;  % Initial point
    x_vals = x(1);  % Store the initial x
    y_vals = x(2);  % Store the initial y
    f_vals = f(x_vals, y_vals);
    grad = grad_f(x(1), x(2));
    k=0;
    while abs(grad)>epsilon
        grad = grad_f(x(1), x(2));
        hess = hessian_f(x(1), x(2));
        %disp(eig(hess))
        eigenvalues = eig(hess);

        isPositiveDefinite = all(eigenvalues > 0);

        if isPositiveDefinite
            dk=-inv(hess)*grad;
            
        else 
            min_e=abs(min(eigenvalues))+1;
            dk=-inv(hess+min_e*[1 0;0 1])*grad;
        end

        %m = max(abs(eigenvalues)) + reg;
        %hess_reg = hess + m*eye(size(hess));
        %disp(x);  % Display the current point
        %d = - gamma * inv(hess_reg) * grad;
        x = x + gamma * dk.';  % Update the point
        % Store the new values of x and y
        x_vals = [x_vals, x(1)];
        y_vals = [y_vals, x(2)];
        f_vals = [f_vals, f(x(1),x(2))];
        
        k=k+1;
        if k>500
            break;
        end
    end
    
    % Final point
    disp('CONSTANT g');
    fprintf("Iterations: %d\n", length(f_vals));
    disp('Final point:');
    disp(['x = ', num2str(x(1))]);
    disp(['y = ', num2str(x(2))]);
    disp(['f(x,y) = ', num2str(f_vals(end))]);
end

function [x_vals, y_vals, f_vals] = lm_method_optimal(f, grad_f, hessian_f, x0, epsilon)
    syms g;
    x = x0;  % Initial point
    x_vals = x(1);  % Store the initial x
    y_vals = x(2);  % Store the initial y
    f_vals = f(x_vals, y_vals);
    grad = grad_f(x(1), x(2));
    k=0;

    while abs(grad)>epsilon
        grad = grad_f(x(1), x(2));
        hess = hessian_f(x(1), x(2));
        %disp(eig(hess))
        eigenvalues = eig(hess);

        isPositiveDefinite = all(eigenvalues > 0);
        if isPositiveDefinite
            dk=-inv(hess)*grad;
        else 
            min_e=abs(min(eigenvalues))+1;
            dk=-inv(hess+min_e*[1 0;0 1])*grad;
        end
        
        %m = max(abs(eigenvalues)) + reg;
        %hess_reg = hess + m*eye(size(hess));
        %disp(x);  % Display the current point

        phi(g) = exp(-(x(1)-g*grad(1))^2 -(x(2)-g*grad(2))^2)*(x(1)-g*grad(1))^5;

        [a,b,~,~] = bisection_method_derivatives(0,8, 0.001, phi);

        % fprintf("#%d: cnt=%d\n",i,cnt);
        
        gamma = (a(end)+b(end))/2;
        %fprintf("Gamma: %d\n", gamma);

        x = x + gamma * dk.';  % Update the point

        % Store the new values of x and y
        x_vals = [x_vals, x(1)];
        y_vals = [y_vals, x(2)];
        f_vals = [f_vals, f(x(1),x(2))];
        k=k+1;
        if k>500
            break;
        end
    end
        
    % Final point
    disp('OPTIMAL g');
    fprintf("Iterations: %d\n", length(f_vals));
    disp('Final point:');
    disp(['x = ', num2str(x(1))]);
    disp(['y = ', num2str(x(2))]);
    disp(['f(x,y) = ', num2str(f_vals(end))]);
end

function [x_vals, y_vals, f_vals] = lm_method_armijo( f, grad_f, hessian_f, x0, epsilon, alpha, beta, g0_armijo)
    x = x0;  % Initial point
    x_vals = x(1);  % Store the initial x
    y_vals = x(2);  % Store the initial y
    f_vals = f(x_vals, y_vals);
    grad = grad_f(x(1), x(2));
    k=0;
    m=0;

    while abs(grad)>epsilon
        grad = grad_f(x(1), x(2));
        hess = hessian_f(x(1), x(2));
        %disp(eig(hess))

        eigenvalues = eig(hess);
        isPositiveDefinite = all(eigenvalues > 0);
        if isPositiveDefinite
            dk=-inv(hess)*grad;
        else 
            min_e=abs(min(eigenvalues))+1;
            dk=-inv(hess+min_e*[1 0;0 1])*grad;
        end


        %m = max(abs(eigenvalues)) + reg;
        %hess_reg = hess + m*eye(size(hess));
        %disp(x);  % Display the current point

        %d = - gamma * inv(hess_reg) * grad;

        % Line search using Armijo rule
        gamma = g0_armijo;
        while f(x(1) + g0_armijo * dk(1)*(beta^m), x(2) + g0_armijo * dk(2)*(beta^m)) > f(x(1), x(2)) + (beta^m)*alpha * g0_armijo * (dk.'* grad) % THE ARMIJO RULE
            m = m + 1;  % Reduce step size
        end
        gamma = g0_armijo * beta^m;

        x = x + gamma * dk.';  % Update the point

        % Store the new values of x and y
        x_vals = [x_vals, x(1)];
        y_vals = [y_vals, x(2)];
        f_vals = [f_vals, f(x(1),x(2))];
        k=k+1;
        if k>500
            break;
        end
    end
    
    % Final point
    disp('ARMIJO g');
    fprintf("Iterations: %d\n", length(f_vals));
    disp('Final point:');
    disp(['x = ', num2str(x(1))]);
    disp(['y = ', num2str(x(2))]);
    disp(['f(x,y) = ', num2str(f_vals(end))]);
end

function [] = tests( f, grad_f, hessian_f, x0, gamma, g0_armijo)
    tol = 1e-6;    % Tolerance for gradient norm
    alpha = 0.001;
    beta = 0.2;
    
    [x_vals_fixed, y_vals_fixed, f_vals_fixed]     = lm_method_fixed(f, grad_f, hessian_f, x0, tol, gamma);
    
    [x_vals_optimal, y_vals_optimal, f_vals_optimal] = lm_method_optimal(f, grad_f, hessian_f, x0, tol);
    
    [x_vals_armijo, y_vals_armijo, f_vals_armijo]   = lm_method_armijo(f, grad_f, hessian_f, x0, tol, alpha, beta, g0_armijo);
    

    % PLOTTING

    % one figure, three plots
    figure;
    tiledlayout(3,1);
    nexttile;
    k_iter = length(f_vals_fixed);
    k = 1:k_iter;
    plot(k, f_vals_fixed, '-o');
    title('CONSTANT g', 'Interpreter', 'latex');
    grid on;
    yline(-0.811174, '--', sprintf('-0.811174'), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    xlabel("Iterations");
    ylabel("f(xk)");
    
    nexttile;
    k_iter = length(f_vals_optimal);
    k = 1:k_iter;
    plot(k, f_vals_optimal, '-o');
    title("OPTIMAL g", 'Interpreter', 'latex');
    grid on;
    yline(-0.811174, '--', sprintf('-0.811174'), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    xlabel("Iterations");
    ylabel("f(xk)");
    
    nexttile;
    k_iter = length(f_vals_armijo);
    k = 1:k_iter;
    plot(k, f_vals_armijo, '-o');
    title("ARMIJO g", 'Interpreter', 'latex');
    grid on;
    yline(-0.811174, '--', sprintf('-0.811174'), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    xlabel("Iterations");
    ylabel("f(xk)");


    % Plot the results in a single 3D plot
    figure;
    
    % Define a grid for plotting the function surface
    [x_grid, y_grid] = meshgrid(linspace(-2, 2, 100), linspace(-2, 2, 100));
    z_grid = f(x_grid, y_grid);  % Evaluate f(x, y) over the grid
    
    % Plot the function surface
    surf(x_grid, y_grid, z_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    colormap parula;
    hold on;
    
    % Specify colors
    color_fixed = [0, 1, 0];  % Green(Fixed Step Size)
    color_optimal = [1, 0, 0];  % Red (Optimal Step Size)
    color_armijo = [0, 0, 1];  % Blue (Armijo Rule)
    
    % Plot the function surface (without adding it to the legend)
    [X, Y] = meshgrid(-2:0.1:2, -2:0.1:2);
    Z = f(X, Y); % Compute function values
    surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Plot the function surface
    hold on;
    
    % Plot each method's trajectory with the specified colors and DisplayName
    h1 = plot3(x_vals_fixed, y_vals_fixed, f(x_vals_fixed, y_vals_fixed), '-o', 'LineWidth', 1.5, 'Color', color_fixed, 'DisplayName', 'Fixed Step Size'); % Fixed step size (Red)
    h2 = plot3(x_vals_optimal, y_vals_optimal, f(x_vals_optimal, y_vals_optimal), '-o', 'LineWidth', 1.5, 'Color', color_optimal, 'DisplayName', 'Optimal Step Size'); % Optimal step size (Green)
    h3 = plot3(x_vals_armijo, y_vals_armijo, f(x_vals_armijo, y_vals_armijo), '-o', 'LineWidth', 1.5, 'Color', color_armijo, 'DisplayName', 'Armijo Rule'); % Armijo rule (Blue)
    
    % Add labels and title
    title(sprintf('[%d, %d]: Trajectories for Different Methods', x0(1), x0(2)));

    xlabel('x');
    ylabel('y');
    zlabel('f(x, y)');
    grid on;
    view(3);  % Set 3D view
    axis tight;  % Adjust axis limits for better visualization
    
    % Add legend to identify methods (exclude the function plot from the legend)
    legend([h1, h2, h3], 'Location', 'northeast', 'FontSize', 12);
end


% 1st TEST
x0 = [0, 0]; % Initial points
gamma = 0.8;
g0_armijo = 8;
tests(f,grad_f,hessian_f, x0, gamma, g0_armijo);
disp(' ');

% 2nd TEST
x0 = [-1, 1]; % Initial points
gamma = 1.5;
g0_armijo = 8;
tests(f,grad_f,hessian_f, x0, gamma, g0_armijo);
disp(' ');

% 3d TEST
x0 = [1, -1]; % Initial points
gamma = 0.5;
g0_armijo = 8;
tests(f,grad_f,hessian_f, x0, gamma, g0_armijo);
disp(' ');
