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
close all;

syms x y z;
f1(x) = (x-2)^2 + x*log(x+3);
f2(y) = exp(-2*y) + (y-2)^2;
f3(z) = exp(z)*(z^3-1) + (z-1)*sin(z);

% 4.1
l = linspace(2*10^(-3), 10^(-1), 20);

count = zeros(length(l), 3);

for i=1:length(l)
    [~,~,~, count(i,1)] = bisection_method_with_derivative(-1, 3, l(i), f1);
    [~,~,~, count(i,2)] = bisection_method_with_derivative(-1, 3, l(i), f2);
    [~,~,~, count(i,3)] = bisection_method_with_derivative(-1, 3, l(i), f3);
end

figure('Name', sprintf('BISECTION METHOD WITH DERIVATIVES A.1'), 'NumberTitle', 'off');
tiledlayout(3,1);
nexttile;
plot(l, count(:,1), '-o');
title('$f_1(x) = (x-2)^2 + xln(x+3)$', 'Interpreter', 'latex');
grid on;
xlabel("l");
ylabel("df1(x)");

nexttile;
plot(l, count(:,2), '-o');
title("$f_2(x) = e^{-2x} + (x-2)^2$", 'Interpreter', 'latex');
grid on;
xlabel("l");
ylabel("df2(x)");

nexttile;
plot(l, count(:,3), '-o');
title("$f_3(x) = e^x(x^3-1) + (x-1)sin(x)$", 'Interpreter', 'latex');
grid on;
xlabel("l");
ylabel("df3(x)");


% 4.2
l = [10^(-3) 10^(-4) 10^(-5) 10^(-6)];

min1 = fminbnd(matlabFunction(f1), -1, 3);
min2 = fminbnd(matlabFunction(f2), -1, 3);
min3 = fminbnd(matlabFunction(f3), -1, 3);

for i=1:length(l)

    [a1, b1, k1, count(i,1)] = bisection_method_with_derivative(-1, 3, l(i), f1);
    [a2, b2, k2, count(i,2)] = bisection_method_with_derivative(-1, 3, l(i), f2);
    [a3, b3, k3, count(i,3)] = bisection_method_with_derivative(-1, 3, l(i), f3);

    figure('Name', sprintf('BISECTION METHOD WITH DERIVATIVES B.%d', i), 'NumberTitle', 'off');
    tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');

    nexttile;
    plot(linspace(1,k1+1,k1+1), a1, '-o');
    hold on;
    grid on;
    plot(linspace(1,k1+1,k1+1), b1, '-o');
    xlabel('Number of Iterations {k}');
    ylabel('ak and bk');
    yline(min1, '--', sprintf('y0 = %.7f', min1), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    legend('ak', 'bk', 'Location', 'best');
    title(sprintf('f1: [ak,bk] VS l = %.1e', l(i)));

    nexttile;
    plot(linspace(1,k2+1,k2+1), a2, '-o');
    hold on;
    grid on;
    plot(linspace(1,k2+1,k2+1), b2, '-o');
    xlabel('Number of Iterations {k}');
    ylabel('ak and bk');
    yline(min2, '--', sprintf('y0 = %.7f', min2), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    legend('ak', 'bk', 'Location', 'best'); 
    title(sprintf('f2: [ak,bk] VS l = %.1e', l(i)));

    nexttile;
    plot(linspace(1,k3+1,k3+1), a3, '-o');
    hold on;
    grid on;
    plot(linspace(1,k3+1,k3+1), b3, '-o');
    xlabel('Number of Iterations {k}');
    ylabel('ak and bk');
    yline(min3, '--', sprintf('y0 = %.7f', min3), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    legend('ak', 'bk', 'Location', 'best'); 
    title(sprintf('f3: [ak,bk] VS l = %.1e', l(i)));
end
