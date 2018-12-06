function xstar = NewtonSolution1D(f,fp,x0,MAX_ITER,EPS)
%
% xstar = NewtonSolution1D(f,fp,x0,MAX_ITER,EPS)
%
% Input:
%    f         :   Function handle of a single variable
%    f         :   Function handle, derivative of f
%    MAX_ITER  :   Maximum number of iteration 
%    EPS       :   Accuracy. (See the description of output)
%    x0        :   Initial guess for the zero of the function f.
%
% Output:
%    xtars     :   Numerical approximation of a zero of f such that
%                  either |f(xstar)| < EPS or the approximation at the
%                  iteration number max_iteration
%
% Last update: February 25, 2018
%
MAX_LUP=1;
xn=x0;
for n = 1:MAX_ITER
    d = -f(xn)/fp(xn); % direction of the solution
    % line search
    Lup = MAX_LUP; % max search distance
    Llo = 0; % minimum search distance
    h = @(lambda) 2*d*f(xn + lambda*d)*fp(xn+lambda*d); 
    % derivative of the minimizer: f^2(xn + lambda*d)
    % looping to minimize h(L)
    for i = 1:MAX_ITER
        L = (Lup + Llo)/2;
        if h(L) > 0
            Lup = L; 
        elseif h(L) < 0
            Llo = L;
        end

        if abs(h(L)) < EPS
            break
        end
    end
    xn = xn + d*L; % approximation after line search
    if abs(f(xn))< EPS
        break
    end
end
xstar = xn;
return
