function u = LinearAdvectionFVSolver1D(U0,a,dx,dt,intervalx,...
    finalT,scheme_option)
%
% u = LinearAdvectionFVSolver1D(U0,a,dx,dt,intervalx,finalT,scheme_option)
%
% Solves u_t + a u_x = 0 with periodic boundary conditions
%
% Input:
%    U0             :   Vectrorized functional handle, antiderivative of 
%                       initial data u(x,0)
%    a              :   Scalar
%    dx             :   The distance between adjacent grid points in space
%    dt             :   The distance between adjacent grid points in time
%    intervalx      :   2 x 1 vector for space domain of u(x,t)
%    finalT         :   A positive number for the final time
%    scheme_option  :   (3) Third order Finite Volume and Runge Kutta
%                           d/dt \bar{u}_j
%                            + 1/Deltax(^f_{j+1/2} - ^f_{j-1/2}) = 0
%                           ^f_{j+1/2} = u^-_{j+1/2} upwind flux
%                           
% Output:
%    u              :   Nt x Nx matrix of numerical approximation of the
%                       solution of Burgers' Equation in 1D. Here Nt and Nx
%                       are the dimensions of discretization space of
%                       intervalx and [0, finalT]
%
% Last update: April 24, 2018

% setting up the discrete problem
xx = intervalx(1) : dx : intervalx(2)-dx;
tt = dt : dt : finalT;
% finite volume scheme
flux1 = @(u) (-1/6)*[u(end) u(1:end-1)] + (5/6)*u + (1/3)*[u(2:end) u(1)];
flux2 = @(u) flux1([u(end) u(1:end-1)]);
f = @(u) (-a/dx)*(flux1(u) - flux2(u)); %RHS of Runge Kutta
Nx = length(xx);
Nt = length(tt);
u = zeros(Nt,Nx);
% initial condition
un = (U0([xx(2:end) xx(1)]) - U0(xx))/dx; % 1 x (Nx - 1)
for nt = 1 : Nt
    switch scheme_option
        case 3 % Runge Kutta 3rd order
            un = RungeKuttaSolver(un,f,nt,dt,3);
    end
    u(nt,:) = un;
end