h = 1/6;
n = 7;
x = 0:h:1;
y = 0:h:1;
z = 0:h:1;

T = zeros(n,n,n);
f_val = -10;

tol = 1e-6;
max_iter = 10000;

for iter = 1:max_iter
    max_diff = 0;
    for i = 2:n-1
        for j = 2:n-1
            for k = 2:n-1
                old = T(i,j,k);
                T(i,j,k) = ( ...
                    T(i+1,j,k) + T(i-1,j,k) + ...
                    T(i,j+1,k) + T(i,j-1,k) + ...
                    T(i,j,k+1) + T(i,j,k-1) ...
                    - h^2 * f_val ) / 6;
                max_diff = max(max_diff, abs(T(i,j,k) - old));
            end
        end
    end
    if max_diff < tol
        break;
    end
end

center = T(4,4,4);
fprintf('T(0.5,0.5,0.5) ≈ %.6f\n', center);

[X,Y] = meshgrid(x,y);

figure;
surf(X,Y,T(:,:,4));
title('Temperature at z = 0.5');
xlabel('x'); ylabel('y'); zlabel('T');