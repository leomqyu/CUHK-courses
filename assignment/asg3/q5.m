k_const = 0.01;
a = 70;
b = 50;
t0 = 0;
tf = 20;
h = 0.5;
t = t0:h:tf;
n = length(t);
y = zeros(1,n);
y(1) = 0;

f = @(t,y) k_const*(a - y)*(b - y);

for i = 1:n-1
    k1 = h * f(t(i), y(i));
    k2 = h * f(t(i) + h/2, y(i) + k1/2);
    k3 = h * f(t(i) + h/2, y(i) + k2/2);
    k4 = h * f(t(i) + h, y(i) + k3);

    y(i+1) = y(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
end

y_exact = 350*(1 - exp(-0.2*t)) ./ (7 - 5*exp(-0.2*t));

error = abs(y_exact - y);

fprintf('%6s %12s %12s %14s\n', 't', 'RK4 y', 'Exact y', 'Abs Error');
fprintf('%s\n', repmat('-',1,48));

for i = 1:n
    fprintf('%6.1f %12.6f %12.6f %14.2e\n', t(i), y(i), y_exact(i), error(i));
end

figure;
plot(t, y, 'b-', 'LineWidth', 1.5); hold on;
plot(t, y_exact, 'r--', 'LineWidth', 1.5);
legend('RK4', 'Exact');
xlabel('t');
ylabel('y(t)');
title('Problem 5: RK4 vs Exact');
grid on;