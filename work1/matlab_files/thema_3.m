% will calculate the first `n` fibonacci numbers of the sequence
% will output a list with these numbers
function fib_sequence = fibonacci(n)
    if n == 1
        fib_sequence = 1;
    elseif n == 2
        fib_sequence = [1, 1];
    else
        fib_sequence = zeros(1, n);
        fib_sequence(1:2) = [1, 1];
        for i = 3:n
            fib_sequence(i) = fib_sequence(i-1) + fib_sequence(i-2);
        end
    end
end

function [alpha, beta, k, cnt] = fibonacci_method(A,B, l, e, f, fibonacci)
    % alpha is ana array that stores a1,a2,...,ai
    % beta the same for b1,b2,...,bi
    % k is the total number of iterations required to finish the algo
    % cnt is the amount of times fi(x) was called to be calculated
    
    alpha = [];
    beta = [];
    n = length(fibonacci(fibonacci <= ((B-A)/l))) + 1; % the +1 includes the Fn>(b-a)/l fibo number
    fibo=fibonacci(1:n);
    x = zeros(2,n);

    if 2*e >= l
        fprintf("The fibonacci method might not work\n");
        cnt = -1;
    else
        alpha(end+1) = A;
        beta(end+1) = B;
        k = 1;
        x(1,k) = alpha(k) + (fibo(n-2)/fibo(n))*(beta(k)-alpha(k));
        x(2,k) = alpha(k) + (fibo(n-1)/fibo(n))*(beta(k)-alpha(k));  
        f1 = f(x(1,k));
        f2 = f(x(2,k));
        cnt = 2;
        while true %breaks through k = n - 2, no need to check |bk-ak|<=l
            if(f1 >= f2)
                %bhma 2
                alpha(end+1) = x(1,k);
                beta(end+1) = beta(k);

                x(1,k+1) = x(2,k); % f(x(2,k)) = f2
                x(2,k+1) = alpha(k+1) + (fibo(n-k-1+1)/fibo(n-k+1))*(beta(k+1)-alpha(k+1)); % every fibo(+1) because of 1-indexed matlab code (F0=fibo(1))
                if k == n-2
                    %bhma 5
                    x(1,n) = x(1,n-1); 
                    x(2,n) = x(1,n-1) + e;
                    cnt = cnt + 1;
                    if f2 >= f(x(2,n)) % f(x(1,n)) = f(x(1,n-1)) = f(x(1,k+1)) = f(x(2,k)) = f2
                        alpha(end+1) = x(1,n);
                        beta(end+1) = beta(n-1);
                        break;
                    else
                        alpha(n) = alpha(n-1);
                        beta(n) = x(2, n);
                        break;
                    end
                else
                    %bhma 4
                    f1 = f2;
                    f2 = f(x(2,k+1));
                    cnt = cnt + 1;
                    k = k+1;
                end
            else 
                % bhma 3
                alpha(end+1) = alpha(k);
                beta(end+1) = x(2, k);
                x(2, k+1) = x(1, k);
                x(1, k+1) = alpha(k+1) + (fibo(n-k-2+1)/fibo(n-k+1))*(beta(k+1)-alpha(k+1)); % every fibo(+1) because of 1-indexed matlab code (F0=fibo(1))
                if k == n-2
                    %bhma 5
                    x(1,n) = x(1,n-1);
                    x(2,n) = x(1,n-1) + e;
                    cnt = cnt + 1;
                    if f1>=f(x(2,n)) % f(x1,n) = f(x1,n-1) = f(x1,k+1) =*testing*= f(x2,k+1) = f(x1,k) = f1
                        alpha(end+1) = x(1,n);
                        beta(end+1) = beta(n-1);
                        break;
                    else
                        alpha(end+1) = alpha(n-1);
                        beta(end+1) = x(2, n);
                        break;
                    end
                else
                    % bhma 4
                    f2 = f1;
                    f1 = f(x(1,k+1));
                    cnt = cnt + 1;
                    k = k+1;
                end
            end
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

% pre-calculate the first 50 fibos and pre-select the Fn you need
% MIGHT CAUSE PROBLEMS FOR very low L -> very high Fn (n>50)
fibo = fibonacci(50);

% 3.1
l = linspace(2*10^(-3), 10^(-1), 20);
count = zeros(length(l), 3);

for i=1:length(l)
    [~,~,~, count(i,1)] = fibonacci_method(-1, 3, l(i), 2*10^(-7), f1, fibo);
    [~,~,~, count(i,2)] = fibonacci_method(-1, 3, l(i), 2*10^(-7), f2, fibo);
    [~,~,~, count(i,3)] = fibonacci_method(-1, 3, l(i), 2*10^(-7), f3, fibo);
end

figure('Name', sprintf('FIBONACCI METHOD A.1'), 'NumberTitle', 'off');
tiledlayout(3,1);
nexttile;
plot(l, count(:,1), '-o');
grid on;
title('$f_1(x) = (x-2)^2 + xln(x+3)$', 'Interpreter', 'latex');
xlabel("l");
ylabel("f1(x)");

nexttile;
plot(l, count(:,2), '-o');
grid on;
title("$f_2(x) = e^{-2x} + (x-2)^2$", 'Interpreter', 'latex');
xlabel("l");
ylabel("f2(x)");

nexttile;
plot(l, count(:,3), '-o');
grid on;
title("$f_3(x) = e^x(x^3-1) + (x-1)sin(x)$", 'Interpreter', 'latex');
xlabel("l");
ylabel("f3(x)");

% 3.2
l = [10^(-3) 10^(-4) 10^(-5) 10^(-6)];

min1 = fminbnd(matlabFunction(f1), -1, 3);
min2 = fminbnd(matlabFunction(f2), -1, 3);
min3 = fminbnd(matlabFunction(f3), -1, 3);

for i=1:length(l)

    [a1, b1, k1, count(i,1)] = fibonacci_method(-1, 3, l(i), 10^(-7), f1, fibo);
    [a2, b2, k2, count(i,2)] = fibonacci_method(-1, 3, l(i), 10^(-7), f2, fibo);
    [a3, b3, k3, count(i,3)] = fibonacci_method(-1, 3, l(i), 10^(-7), f3, fibo);

    figure('Name', sprintf('FIBONACCI METHOD B.%d', i), 'NumberTitle', 'off');
    tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');

    nexttile;
    plot(linspace(1,k1+2,k1+2), a1, '-o'); % EINAI EIDIKO GIA TON FIBO ALGORITHMO
    hold on;
    grid on;
    plot(linspace(1,k1+2,k1+2), b1, '-o');
    xlabel('Number of Iterations {k}');
    ylabel('ak and bk');
    yline(min1, '--', sprintf('y0 = %.7f', min1), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    legend('ak', 'bk', 'Location', 'best');
    title(sprintf('f1: [ak,bk] VS l = %.1e', l(i)));

    nexttile;
    plot(linspace(1,k2+2,k2+2), a2, '-o');
    hold on;
    grid on;
    plot(linspace(1,k2+2,k2+2), b2, '-o');
    xlabel('Number of Iterations {k}');
    ylabel('ak and bk');
    yline(min2, '--', sprintf('y0 = %.7f', min2), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    legend('ak', 'bk', 'Location', 'best'); 
    title(sprintf('f2: [ak,bk] VS l = %.1e', l(i)));

    nexttile;
    plot(linspace(1,k3+2,k3+2), a3, '-o');
    hold on;
    grid on;
    plot(linspace(1,k3+2,k3+2), b3, '-o');
    xlabel('Number of Iterations {k}');
    ylabel('ak and bk');
    yline(min3, '--', sprintf('y0 = %.7f', min3), 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
    legend('ak', 'bk', 'Location', 'best'); 
    title(sprintf('f3: [ak,bk] VS l = %.1e', l(i)));
end
