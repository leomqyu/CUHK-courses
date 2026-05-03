tol = 1e-6;
max_iter = 50;
x = 1;

for iter = 1:max_iter
    f_integrand = @(t) exp(x * t.^2);
    fp_integrand = @(t) t.^2 .* exp(x * t.^2);
    fx = adaptive_simpson(f_integrand, 0, 1, 1e-6) - 2;
    fpx = adaptive_simpson(fp_integrand, 0, 1, 1e-6);
    fprintf('At iteration %d:\n', iter);
    fprintf('x = %.6f\n', x);
    fprintf('f = %.6f\n', fx);
    fprintf("f' = %.6f\n\n", fpx);
    x_new = x - fx / fpx;
    if abs(x_new - x) < tol
        break;
    end
    x = x_new;
end

fprintf('Root by Newton Method: %.6f\n', x);

function I = adaptive_simpson(f, a, b, tol)
    mid = (a + b)/2;
    
    I1 = simpson(f, a, b);
    I2 = simpson(f, a, mid) + simpson(f, mid, b);
    
    if abs(I2 - I1) < tol
        I = I2;
    else
        I = adaptive_simpson(f, a, mid, tol/2) + adaptive_simpson(f, mid, b, tol/2);
    end
end

function I = simpson(f, a, b)
    h = (b - a)/2;
    I = (h/3) * (f(a) + 4*f(a + h) + f(b));
end