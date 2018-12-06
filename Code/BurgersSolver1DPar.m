function u = BurgersSolver1DPar(u0,u0p,dx,dt,intervalx,intervalt,MAX_ITER,EPS)
%
% u = BurgersSolver1DPar(u0,u0p,dx,dt,intervalx,intervalt,MAX_ITER,EPS)
%
% Input:
%    u0        :   Functional handle for initial data: u(x,0)
%    u0p       :   Function handle, derivative of u0
%    dx        :   The distance between adjacent grid points in space                
%    dt        :   The distance between adjacent grid points in time
%    intervalx :   2 x 1 vector for space domain of u(x,t)
%    intervalt :   2 x 1 vector for time domain of u(x,t)
%    MAX_ITER  :   An integer. Maximum number of iterations that is allowed
%                  for the Newton's iteration (Suggested 20)
%    EPS       :   Accuracy of Newton's method. (Suggested 1e-05)
%
% Output:
%    u         :   Nt x Nx matrix of numerical approximation of the 
%                  solution of Burgers' Equation in 1D. Here Nt and Nx are
%                  the dimensions of discretization space of inervalx and
%                  intervaly
%
% Last update: March 13, 2018

% setting up the discrete problem
xx = intervalx(1) : dx : intervalx(2);
tt = intervalt(1)+dt : dt : intervalt(2);
% looping over x&t to find calculate solution u(x,t)
Nx = length(xx);
Nt = length(tt);
u = zeros(Nt,Nx);
xx1=xx(1);
parfor nt = 1 : Nt
    v = zeros(1,Nx);
    x0 = xx1;
    f = @(x) tt(nt) * u0 (x) + x - xx1;
    fp = @(x) tt(nt) * u0p (x) + 1;
    % using Newton's method to find the solution to
    % (x - xstar)/(t - 0) = u0(xstar)
    xstar = NewtonSolution1D(f,fp,x0,MAX_ITER,EPS);
    v(1) = u0(xstar);
    for nx = 2 : Nx
        f = @(x) tt(nt) * u0 (x) + x - xx(nx);
        fp = @(x) tt(nt) * u0p (x) + 1;
        x0 = xstar;
        % using Newton's method to find the solution to
        % (x - xstar)/(t - 0) = u0(xstar)
        xstar = NewtonSolution1D(f,fp,x0,MAX_ITER,EPS);
        v(nx) = u0(xstar);
    end
    u(nt,:) = v;
end