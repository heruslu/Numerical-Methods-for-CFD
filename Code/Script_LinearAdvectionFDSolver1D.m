% Script to test LinearAdvection1DSolver
% Last update: March 13, 2018
% Solves u_t + a u_x = 0 with periodic boundary conditions
% initial condition
u0 = @(x) sin(2*pi*x);
a=1;
display_option = 1;%(1) error for given spatial refinement
%(2) plotting an animated solution
scheme_option = 1;%(1) u^(n+1)_(j) = u^_j - a*c*(u^_j - u^_(j-1))%notstable
%(2) u^(n+1)_(j) = u^_j - a*c*(u^_(j+1) - u^_j)% stable
% discretization parameteres
intervalx = [0 1];
%dx_vector=1./[10,20,40,80,160,320,640,1280]; % space discretization
dx_vector=1./[20,40,80,160,320,640,1280]; % space discretization
finalT = 3; % the final time
CFL=1;
dt_vector=CFL*dx_vector; % time discretization
% solver
err_vector=[];
for ind=1:length(dx_vector)
    dx=dx_vector(ind);
    dt=dt_vector(ind);
    u = LinearAdvectionSolver1D(u0,a,dx,dt,intervalx,finalT,...
        scheme_option);
    xx = intervalx(1) : dx : intervalx(2);
    uexac=u0(xx-a*finalT); % exact solution at final time
    error=dx*sum(abs(u(end,:)-uexac));
    err_vector=[err_vector;error];%#ok
end
%
switch display_option
    case 1
        disp('Error and the order of convergence:')
        order = log(err_vector(1:end-1)./err_vector(2:end))/log(2);
        order = [0;order];
        disp([err_vector order]);
    case 2
        disp('Solution for the refined mesh')
        tt = 0+dt : dt : finalT-dt;
        figure
        for nt = 1 : size(u,1)
            plot(xx, u(nt,:),xx, u0(xx-a*tt(nt)));
            title(['time = ' num2str(tt(nt))]);
            pause(0.01);
        end
end