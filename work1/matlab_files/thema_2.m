function [alpha, beta, k, cnt] = golden_section_method(A,B, l, f)
    % alpha is ana array that stores a1,a2,...,ai
    % beta the same for b1,b2,...,bi
    % k is the total number of iterations required to finish the algo
    % cnt is the amount of times fi(x) was called to be calculated

    gama = (-1+sqrt(5))/2; % phi calculated
    alpha = [];
    beta = [];
    k = 1;

    alpha(end+1) = A;
    beta(end+1) = B;
    % pre-calculations
    x1 = alpha(k) + (1-gama)*(beta(k)-alpha(k));
    x2 = alpha(k) + gama*(beta(k)-alpha(k));    
    f1 = f(x1);
    f2 = f(x2);
    cnt = 2;
    while abs(beta(k)-alpha(k))>=l
        if(f1 >= f2)
            alpha(end+1) = x1;
            beta(end+1) = beta(k);

            x1 = x2;
            f1 = f2; % SKIPPING ONE F(X) CALCULATION

            x2 = alpha(k+1) + gama*(beta(k+1)-alpha(k+1));
            f2 = f(x2);
        else 
            alpha(end+1) = alpha(k);
            beta(end+1) = x2;

            x2 = x1;
            f2 = f1; % SKIPPING ONE F(X) CALCULATION

            x1 = alpha(k+1) + (1-gama)*(beta(k+1)-alpha(k+1));
            f1 = f(x1);
        end
        cnt = cnt + 1;
        k = k+1;
    end
end

clearvars; % clear workspace
clc; % clear command window
close all;

syms x y z;
f1(x) = (x-2)^2 + x*log(x+3);
f2(y) = exp(-2*y) + (y-2)^2;
f3(z) = exp(z)*(z^3-1) + (z-1)*sin(z);

% 2.1
l = linspace(2*10^(-3), 10^(-1), 20);

count = zeros(length(l), 3);

for i=1:length(l)
    [~,~,~, count(i,1)] = golden_section_method(-1, 3, l(i), f1);
    [~,~,~, count(i,2)] = golden_section_method(-1, 3, l(i), f2);
    [~,~,~, count(i,3)] = golden_section_method(-1, 3, l(i), f3);
end

figure('Name', sprintf('GOLDEN SECTION METHOD A.1'), 'NumberTitle', 'off');
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


% 2.2
l = [10^(-3) 10^(-4) 10^(-5) 10^(-6)];

min1 = fminbnd(matlabFunction(f1), -1, 3);
min2 = fminbnd(matlabFunction(f2), -1, 3);
min3 = fminbnd(matlabFunction(f3), -1, 3);

for i=1:length(l)

    [a1, b1, k1, ~] = golden_section_method(-1, 3, l(i), f1);
    [a2, b2, k2, ~] = golden_section_method(-1, 3, l(i), f2);
    [a3, b3, k3, ~] = golden_section_method(-1, 3, l(i), f3);

    figure('Name', sprintf('GOLDEN SECTION METHOD B.%d', i), 'NumberTitle', 'off');
    tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');

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
