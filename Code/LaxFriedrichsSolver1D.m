function uh = LaxFriedrichsSolver1D(U0,f,fp,dx,dt,intervalx,T,...
    scheme_option)
%
% uh = LaxFriedrichsSolver1D(U0,f,fp,dx,dt,intervalx,T,scheme_op)
%
% Solves u_t + f(u)_x = 0    for x in intervalx, t in [0,T]
% with periodic boundary conditions
%
% Input:
%    U0            :   Function handle, antiderivative of initial data u0
%    f             :   Function handle
%    fp            :   Function handle, derivative of f
%    dx            :   The distance between adjacent grid points in space
%    dt            :   The distance between adjacent grid points in time
%    intervalx     :   2 x 1 vector for space domain of u(x,t)
%    T             :   Final time for time domain of u(x,t)
%    scheme_option :   A number from the set {1,3,4} where
%                      (1) 1st space and time
%                      (3) 3rd space and time scheme
%                      (4) 3rd space and time scheme with minmod correction
%
% Output:
%    uh            :   Nt x Nx matrix of numerical approximation of the
%                      solution of the PDE in 1D. Here Nt and Nx are the
%                      dimensions of discretization spaces for t and x
%
% Last update: December 4, 2018

% setting up the discrete problem
xx = intervalx(1) : dx : intervalx(2)-dx;
tt = dt : dt : T;
% initial condition
un = (U0([xx(2:end) xx(1)]) - U0(xx))/dx; % 1 x (Nx - 1)
% initializing the solution
Nx = length(xx); Nt = length(tt);
uh = zeros(Nt,Nx);
alpha = max(abs(fp(un)));
uminus = @(u) (-1/6)*[u(end) u(1:end-1)] + (5/6)*u + (1/3)*[u(2:end) u(1)];
uplus = @(u) (1/3)*u + (5/6)*[u(2:end) u(1)] - (1/6)*[u(3:end) u(1:2)];
switch scheme_option
    case 1 % 1st order space and time
        flux = @(u)0.5*(f(u)+f([u(2:end) u(1)])+alpha*(u-[u(2:end) u(1)]));
    case 3 % 3rd space and time
        flux = @(u) 0.5*(f(uminus(u))+f(uplus(u))+...
            alpha*(uminus(u)-uplus(u)));
    case 4 % 3rd order space and time with minmod correction
        signer = @(a,b,c) 1-0.5*(abs(sign(a)-sign(b)) + ...
            abs(sign(b)-sign(c)) + abs(sign(c)-sign(a)));
        minmod = @(a,b,c) sign(a).*signer(a,b,c)...
            .*min([abs(a);abs(b);abs(c)],[],1);
        uminusUp = @(u) u + minmod(uminus(u)-u,[u(2:end) u(1)]-u,...
            u-[u(end) u(1:end-1)]);
        uplusUp = @(u) [u(2:end) u(1)] - minmod([u(2:end) u(1)] - ...
            uplus(u),[u(2:end) u(1)]-u,[u(3:end) u(1:2)]-[u(2:end) u(1)]);
        flux = @(u) 0.5*(f(uminusUp(u))+f(uplusUp(u))+...
            alpha*(uminusUp(u)-uplusUp(u)));
end

spatial_disc = @(yn) (-1/dx)*(flux(yn) - flux([yn(end) yn(1:end-1)]));
linear_scheme=(scheme_option==1);
for nt = 1 : Nt
    if ~linear_scheme
        un = RungeKuttaSolver(un,spatial_disc,dt,3);
    else
        un = un + dt*spatial_disc(un);
    end
    uh(nt,:) = un;
end