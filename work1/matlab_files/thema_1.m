function [alpha, beta, k, cnt] = bisection_method(A,B, l, e, f)
    % alpha is ana array that stores a1,a2,...,ai
    % beta the same for b1,b2,...,bi
    % k is the total number of iterations required to finish the algo
    % cnt is the amount of times fi(x) was called to be calculated

    alpha = [];
    beta = [];
    k = 1;
    cnt = 0;

    alpha(end+1) = A;
    beta(end+1) = B;
    while abs(beta(k)-alpha(k))>l
        x1 = (alpha(k)+beta(k))/2 - e;
        x2 = (alpha(k)+beta(k))/2 + e;

        if(f(x1)>=f(x2)) % 2 calculations of f(x) each iteration
            alpha(end+1) = x1;
            beta(end+1) = beta(k);
        else 
            alpha(end+1) = alpha(k);
            beta(end+1) = x2;
        end
        cnt = cnt + 2;
        k = k+1;
    end
end

clearvars; % clear workspace
clc;       % clear command window
close all;

syms x y z;
f1(x) = (x-2)^2 + x*log(x+3);
f2(y) = exp(-2*y) + (y-2)^2;
f3(z) = exp(z)*(z^3-1) + (z-1)*sin(z);

% 1.1
e = linspace(10^(-5), 10^(-2)/2.1, 20); % specified number of points 
%e = linspace(10^(-6), 0.0099/2, 20); % specified number of points 
l = 10^(-2);
count = zeros(length(e), 3);

for i=1:length(e)
    [~, ~, ~, count(i,1)] = bisection_method(-1, 3, l, e(i), f1);
    [~, ~, ~, count(i,2)] = bisection_method(-1, 3, l, e(i), f2);
    [~, ~, ~, count(i,3)] = bisection_method(-1, 3, l, e(i), f3);
end

% one figure, three plots and a little bit of sauce
tiledlayout(3,1);
nexttile;
plot(e, count(:,1), '-o');
title('$f_1(x) = (x-2)^2 + xln(x+3)$', 'Interpreter', 'latex');
grid on;
xlabel("e");
ylabel("f1(x)");

nexttile;
plot(e, count(:,2), '-o');
title("$f_2(x) = e^{-2x} + (x-2)^2$", 'Interpreter', 'latex');
grid on;
xlabel("e");
ylabel("f2(x)");

nexttile;
plot(e, count(:,3), '-o');
title("$f_3(x) = e^x(x^3-1) + (x-1)sin(x)$", 'Interpreter', 'latex');
grid on;
xlabel("e");
ylabel("f3(x)");


% 1.2
e = 10^(-3);
l = linspace(2.001*10^(-3), 10^(-1), 20);
count = zeros(length(l), 3);
for i=1:length(l)
    [~,~,~, count(i, 1)] = bisection_method(-1, 3, l(i), e, f1);
    [~,~,~, count(i, 2)] = bisection_method(-1, 3, l(i), e, f2);
    [~,~,~, count(i, 3)] = bisection_method(-1, 3, l(i), e, f3);
end

% one figure, 3 plots, and some sauce
figure;
tiledlayout(3,1);
nexttile;
plot(l, count(:, 1), '-o');
title('$f_1(x) = (x-2)^2 + xln(x+3)$', 'Interpreter', 'latex');
grid on;
xlabel("l");
ylabel("f1(x)");

nexttile;
plot(l, count(:, 2), '-o');
title("$f_2(x) = e^{-2x} + (x-2)^2$", 'Interpreter', 'latex');
grid on;
xlabel("l");
ylabel("f2(x)");

nexttile;
plot(l, count(:, 3), '-o');
title("$f_3(x) = e^x(x^3-1) + (x-1)sin(x)$", 'Interpreter', 'latex');
grid on;
xlabel("l");
ylabel("f3(x)");


% 1.3 
e = 10^(-6)/2.001; % gia na isxuei 2e<lmin
l = [10^(-3) 10^(-4) 10^(-5) 10^(-6)];

min1 = fminbnd(matlabFunction(f1), -1, 3);
min2 = fminbnd(matlabFunction(f2), -1, 3);
min3 = fminbnd(matlabFunction(f3), -1, 3);

for i=1:length(l)
    [a1, b1, k1, ~] = bisection_method(-1, 3, l(i), e, f1);
    [a2, b2, k2, ~] = bisection_method(-1, 3, l(i), e, f2);
    [a3, b3, k3, ~] = bisection_method(-1, 3, l(i), e, f3);
    
    % code to have 3 plots on one figure, and create one figure per
    % iteration of main loop
    figure;
    tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');

    
    % standard code for plotting and 
    % adding a little bit of sauce
    nexttile;
    plot(linspace(1,k1,k1), a1, '-o');
    hold on;
    grid on;
    plot(linspace(1,k1,k1), b1, '-o');
    xlabel('Number of Iterations {k}');
    ylabel('ak and bk');
    yline(min1, '--', sprintf('y0 = %.7f', min1), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    legend('ak', 'bk', 'Location', 'best'); 
    title(sprintf('f1: [ak,bk] VS l = %.1e', l(i)));

    % ctrl+c ctrl+v
    nexttile;
    plot(linspace(1,k2,k2), a2, '-o');
    hold on;
    grid on;
    plot(linspace(1,k2,k2), b2, '-o');
    xlabel('Number of Iterations {k}');
    ylabel('ak and bk');
    yline(min2, '--', sprintf('y0 = %.7f', min2), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    legend('ak', 'bk', 'Location', 'best'); 
    title(sprintf('f2: [ak,bk] VS l = %.1e', l(i)));

    % ctrl+c ctrl+v
    nexttile;
    plot(linspace(1,k3,k3), a3, '-o');
    hold on;
    grid on;
    plot(linspace(1,k3,k3), b3, '-o');
    xlabel('Number of Iterations {k}');
    ylabel('ak and bk');
    yline(min3, '--', sprintf('y0 = %.7f', min3), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    legend('ak', 'bk', 'Location', 'best'); 
    title(sprintf('f3: [ak,bk] VS l = %.1e', l(i)));
end