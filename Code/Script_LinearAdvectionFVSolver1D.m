% Script to test LinearAdvectionFVSolver1D
% Last update: April 24, 2018
% Solves u_t + a u_x = 0 with periodic boundary conditions
% Contains: Third order time and space discretization

% PARAMETERS
a=1; % coefficient of u_x
intervalx = [-pi pi]; % interval of x
finalT = 1; % the final time
% initial condition
u0 = @(x) sin(x);
U0 = @(x) -cos(x); % integral of u0
scheme_option = 3; % 3rd order Runge Kutta (Currenlty only one option)
% discretization parameteres
refinement_size = [40,80,160,320,640,1280,5120];
% space discretization
dx_vector=(intervalx(2)-intervalx(1))./refinement_size; 
CFL=1; dt_vector=CFL*dx_vector; % time discretization
% display setting
display_option = 1;%(1) error for given spatial refinement
                   %(2) plotting an animated solution

% SOLVER
err_vector=[];
for ind=1:length(refinement_size)
    dx=dx_vector(ind);
    dt=dt_vector(ind);
    xx = intervalx(1) : dx : intervalx(2)-dx;
    u = LinearAdvectionFiniteVolumeSolver1D(U0,a,dx,dt,intervalx,finalT,...
        scheme_option);
    % exact solution at final time
    numFinalT = dt*floor(finalT/dt); % approximation to final time
    uexac = (U0([xx(2:end) xx(1)]-numFinalT) ...
        - U0(xx-numFinalT))/dx; % 1 x (Nx - 1)
    error=sum(abs(u(end,:)-uexac))/sum(abs(uexac));
    err_vector=[err_vector;error];%#ok
end

% DISPLAY
disp('Error and the order of convergence:')
order = log(err_vector(1:end-1)./err_vector(2:end))/log(2);
order = [0;order];
disp([err_vector order]);
switch display_option
    case 2
        disp('Solution for the refined mesh')
        tt = dt : dt : finalT;
        figure
        for nt = 1 : size(u,1)
            uexac = (U0([xx(2:end) xx(1)]-numFinalT) ...
                - U0(xx-numFinalT))/dx; % 1 x (Nx - 1)
            plot(xx, u(nt,:), xx, uexac);
            title(['time = ' num2str(tt(nt))]);
            pause(0.01);
        end
        legend('Finite Volume','Exact');
end