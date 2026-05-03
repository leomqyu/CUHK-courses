A = [ 2  3  1; 4  7  8; -2  4  5 ];
b = [1; 2; 3];
x = solveLinSys(A, b);
fprintf('solution x:\n');
fprintf('%g\n', x);

% % Check 
% residual = A*x - b;
% disp(norm(residual));


function x = solveLinSys(A, b)
    % do LU factorization
    [L, U, P] = PartialPivotLU(A);

    % permute b
    b_tilde = P * b;

    % cal y
    y = forwardSubSAXPY(L, b_tilde);

    % cal x
    x = backSubSAXPY(U, y);
end

function [L, U, P] = PartialPivotLU(A)
    [n, m] = size(A);
    if n ~= m, error(''); end

    U = A;
    L = eye(n);
    P = eye(n);

    for k = 1:n-1
        % find pivot for kth col
        [~, p] = max(abs(U(k:n, k)));
        p = p + k - 1;

        if p ~= k
            % swap
            U([k p], :) = U([p k], :);
            P([k p], :) = P([p k], :);
            if k > 1
                L([k p], 1:k-1) = L([p k], 1:k-1);
            end
        end

        % cal L,U
        for i = k+1:n
            L(i, k) = U(i, k) / U(k, k);
            U(i, :) = U(i, :) - L(i, k) * U(k, :);
        end
    end
end

function x = backSubSAXPY(U, b)
    n = length(b);
    for i = n:-1:1
        x(i) = b(i)/U(i,i);
        b(1:i-1) = b(1:i-1) - x(i)*U(1:i-1,i);
    end
    x=x(:);
end

function x = forwardSubSAXPY(L, b)
    n = length(b);
    for i = 1:n
        x(i) = b(i);                 % since L(i,i) = 1
        b(i+1:n) = b(i+1:n) - x(i)*L(i+1:n, i);
    end
    x = x(:);
end