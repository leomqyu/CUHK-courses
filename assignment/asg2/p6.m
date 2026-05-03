% (a)
function evals = myShiftedQR(A, maxIter, tol)
    n = size(A, 1);
    evals = zeros(n, 1);
    Ak = A;
    
    for k = n:-1:2
        for iter = 1:maxIter
            mu = Ak(k, k); 
            [Q, R] = qr(Ak(1:k, 1:k) - mu * eye(k));
            Ak(1:k, 1:k) = R * Q + mu * eye(k);
            if abs(Ak(k, k-1)) < tol
                break;
            end
        end
        evals(k) = Ak(k, k);
    end
    evals(1) = Ak(1, 1);
end

function [lambda, x, errors] = myPowerIteration(A, tol)
    n = size(A, 1);
    x = rand(n, 1);
    x = x / norm(x);
    lambda = 0;
    errors = [];
    maxIter = 1000;
    
    for i = 1:maxIter
        x_new = A * x;
        lambda_new = (x_new' * x) / (x' * x); 
        errors(i) = abs(lambda_new - 1);
        if abs(lambda_new - lambda) < tol
            break;
        end
        x = x_new / norm(x_new);
        lambda = lambda_new;
    end
end

fprintf('\n(a)\n');
rng(42);
A = rand(20);
A = A ./ sum(A); 

matlab_evals = eig(A);
[~, idx] = max(abs(matlab_evals)); 
lambda_true = matlab_evals(idx); % The dominant eigenvalue
qr_evals = myShiftedQR(A, 1000, 1e-12);
[~, qr_idx] = max(abs(qr_evals));
lambda_qr = qr_evals(qr_idx);
[lambda_pow, ~, ~] = myPowerIteration(A, 1e-12);

diff_qr = abs(lambda_qr - lambda_true);
diff_pow = abs(lambda_pow - lambda_true);

fprintf('MATLAB eig() dominant eigenvalue:  %.15f\n', lambda_true);
fprintf('myShiftedQR dominant eigenvalue:   %.15f\n', lambda_qr);
fprintf('myPowerIteration dominant eigenvalue: %.15f\n', lambda_pow);

% (b)
tol = 1e-10;
maxIter = 1000;
[~, ~, err_pow] = myPowerIteration(A, tol);
err_qr = [];
Ak = A;
for iter = 1:maxIter
    mu = Ak(end, end); 
    [Q, R] = qr(Ak - mu * eye(size(Ak)));
    Ak = R * Q + mu * eye(size(Ak));
    current_evals = diag(Ak);
    [~, idx] = min(abs(current_evals - 1));
    err_qr(iter) = abs(current_evals(idx) - 1);
    if err_qr(iter) < tol, break; end
end
figure;
semilogy(err_pow, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
semilogy(err_qr, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 4);
grid on;
xlabel('Iteration Number');
ylabel('Error');
title('Error against Iteration Number');
legend('Power Iteration', 'Shifted QR');