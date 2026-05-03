T = 500000;
N_vals = zeros(T,1);

for t = 1:T
    sum_val = 0;
    count = 0;
    while sum_val <= 1
        sum_val = sum_val + rand();
        count = count + 1;
    end
    N_vals(t) = count;
end

mu_hat = mean(N_vals);
sigma_hat = std(N_vals);

true_val = exp(1);
abs_error = abs(mu_hat - true_val);

CI_low = mu_hat - 1.96 * sigma_hat / sqrt(T);
CI_high = mu_hat + 1.96 * sigma_hat / sqrt(T);

fprintf('(1) Estimate: %.5f\n', mu_hat);
fprintf('(2) Absolute Error: %.5f\n', abs_error);
fprintf('(3) 95%% CI: [%.5f, %.5f]\n', CI_low, CI_high);