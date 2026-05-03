% (a)
function c_dft = dft_direct(y)
    n = length(y);
    c_dft = zeros(n, 1);
    for k = 0:n-1
        sum_val = 0;
        for j = 0:n-1
            sum_val = sum_val + y(j+1) * exp(-2 * pi * 1i * j * k / n);
        end
        c_dft(k+1) = sum_val / n;
    end
end

% (b)
function p = triginterp_eval(c, xq)
    n = length(c);
    p = zeros(size(xq));
    
    ks = (-n/2):(n/2);
    
    c_shifted = [c(n/2+1); c(n/2+2:end); c(1); c(2:n/2); c(n/2+1)];
    
    weights = ones(size(ks));
    weights(1) = 0.5;
    weights(end) = 0.5;
    
    for i = 1:length(xq)
        p(i) = sum(weights' .* c_shifted .* exp(1i * ks' * xq(i)));
    end
end


% (c)
n_c = 16;
f = @(x) exp(sin(x));
xj = 2 * pi * (0:n_c-1)' / n_c;
yj = f(xj);

c = dft_direct(yj);
p_vals = triginterp_eval(c, xj);

max_err = max(abs(p_vals - yj));
fprintf('(c): Max error for n=16: %.4e\n', max_err);

% (d)
rng(42);
m_vals = [4, 6, 8, 10, 12];
n_vals = 2.^m_vals;

fprintf('\n(d)\n');
fprintf('%-6s | %-15s | %-12s | %-12s\n', 'n', 'Rel. Inf Error', 'Direct Time', 'FFT Time');
fprintf('-------------------------------------------------------------\n');

for n = n_vals
    y_rand = rand(n, 1);
    
    tic;
    c_direct = dft_direct(y_rand);
    t_direct = toc;
    
    tic;
    c_fft = fft(y_rand) / n; 
    t_fft = toc;
    
    rel_error = norm(c_direct - c_fft, Inf) / norm(c_fft, Inf);
    
    fprintf('%-6d | %-15.4e | %-12.4f | %-12.5f\n', n, rel_error, t_direct, t_fft);
end