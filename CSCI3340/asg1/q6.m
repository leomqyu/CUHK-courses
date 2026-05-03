% init
x0 = 1;           
tol = 1e-10;  

f = @(x) 4*x^3 - 52*x^2 + 160*x - 100;
% f = @(x) (16-2*x)*(10-2*x)*x - 100;

% Newton's method
x_k = x0;

while true
    f_xk = f(x_k);
    
    J = numJacobian(f, x_k);

    x_next = x_k - f_xk/J;   % for scalor J
    
    if abs(x_next - x_k) < tol
        s = x_next;
        break;
    end
    
    x_k = x_next;
end

fprintf('Solution of f(x)=0: x=%.5f\n', s)


function J = numJacobian(f, x)
    for i= 1:length(x)
        xd = x;
        h = sqrt(eps*(1+abs(xd(i))));
        xd(i) = xd(i) + h;
        J(:, i) = (f(xd)-f(x))/h;
    end
end