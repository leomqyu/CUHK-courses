w = [48; 61; 81; 113; 131];
p = [91; 98; 103; 110; 112];
A = [ones(size(w)), log(w)];
b = p;

% (b)
fprintf('\n(b)\n');
[U, S, V] = svd(A, "econ");
S_inv = diag(1 ./ diag(S));
A_plus = V * S_inv * U';
beta_hat = A_plus * b;

sigma = diag(S);
fprintf('Beta_0: %.4f, Beta_1: %.4f\n', beta_hat(1), beta_hat(2));
fprintf('Singular values: %.4f, %.4f\n', sigma(1), sigma(2));

% (c)
fprintf('\n(c)\n');
r = b - A * beta_hat;
norm_r = norm(r, 2);
rel_norm_r = norm_r / norm(b, 2);
fprintf('norm_r: %.4f, rel_norm_r: %.4f\n', norm_r, rel_norm_r);

w_fine = linspace(40, 150, 100);
p_fit = beta_hat(1) + beta_hat(2) * log(w_fine);
figure;
plot(w, p, 'ro', w_fine, p_fit, 'b-');
grid on; legend('Data Points', 'Fitted Curve');

% (d) 
fprintf('\n(d)\n');
kappa_A = sigma(1) / sigma(2);
b_tilde = b + 0.1 * randn(size(b));
beta_tilde = A_plus * b_tilde;
rel_beta_error = norm(beta_tilde - beta_hat, 2) / norm(beta_hat, 2);
fprintf('Condition Number: %.4f\n', kappa_A);
fprintf('Relative error in Beta: %.4e\n', rel_beta_error);