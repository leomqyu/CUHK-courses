x = 3;
fprintf('f(%d) = %.6f\n', x, f(x));

function y = f(x, tol)
    % find tol
    if(nargin < 2)
        tol = 0;
    end

    % calculate until stop
    n = 1;
    t_n = x^2;
    y = t_n;

    while(t_n > tol)
        n = n + 1;
        t_n = (2*n^2-1)/((2*(n-1)^2-1)*n) * x^2 * t_n;
        y = y + t_n;
    end
end