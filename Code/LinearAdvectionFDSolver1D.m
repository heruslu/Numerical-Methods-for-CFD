function u = LinearAdvectionFDSolver1D(u0,a,dx,dt,intervalx,finalT,...
    scheme_option)
%
% u = LinearAdvectionFDSolver1D(u0,a,dx,dt,intervalx,finalT,scheme_option)
%
% Solves u_t + c u_x = 0 with periodic boundary conditions
%
% Input:
%    u0             :   Vectrorized functional handle for initial data
%                       u(x,0)
%    a              :   Scalar
%    dx             :   The distance between adjacent grid points in space
%    dt             :   The distance between adjacent grid points in time
%    intervalx      :   2 x 1 vector for space domain of u(x,t)
%    finalT         :   A positive number for the final time
%    scheme_option  :   (1) u^(n+1)_(j) = u^_j - c*lambda*(u^_j - u^_(j-1))
%                       (2) u^(n+1)_(j) = u^_j - c*lambda*(u^_(j+1) - u^_j)
% Output:
%    u              :   Nt x Nx matrix of numerical approximation of the
%                       solution of Burgers' Equation in 1D. Here Nt and Nx
%                       are the dimensions of discretization space of
%                       intervalx and [0, finalT]
%
% Last update: March 13, 2018

intervalt = [0 finalT];
% setting up the discrete problem
xx = intervalx(1) : dx : intervalx(2);
tt = intervalt(1)+dt : dt : intervalt(2)-dt;
% looping over x&t to find calculate solution u(x,t)
Nx = length(xx);
Nt = length(tt);
u = zeros(Nt,Nx);
un = u0(xx);
c = dt/dx;
for nt = 1 : Nt
    switch scheme_option
        case 1 % (1) u^(n+1)_(j) = u^n_j - c*lambda*(u^n_j - u^n_(j-1))
            un(2:end) = un(2:end) - a*c*(un(2:end) - un(1:end-1)); % 1 x Nx vector
            un(1) = un(end);    
        case 2 % (2) u^(n+1)_(j) = u^n_j - c*lambda*sw(u^n_(j+1) - u^n_j)
            un(1:end-1) = un(1:end-1) - a*c*(un(2:end) - un(1:end-1));
            un(end) = un(1);
    end
    u(nt,:) = un;
end