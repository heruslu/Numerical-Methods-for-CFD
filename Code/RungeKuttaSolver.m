function ynp1 = RungeKuttaSolver(yn,f,dt,order)
%
% ynp1 = RungeKuttaSolver(yn,f,dt,order)
%
% Implements one step of Runge Kutta method for the ODE
% d/dt y = f(y,t)
%
%
% Input:
%    yn     :   N x 1 vector, current time step approximation
%    f      :   Vectorized function handle
%    dt     :   Real number, step length
%    order  :   Integer
%
% Output:
%    ynp1   :   N x 1 vector, next time step approximation
%
% Last update: April 23, 2018

switch order
    case 3
        y1 = yn + dt*f(yn);
        y2 = (3/4)*yn + (1/4)*y1 + (dt/4)*f(y1);
        ynp1 = (1/3)*yn + (2/3)*y2 + (2*dt/3)*f(y2);
end
