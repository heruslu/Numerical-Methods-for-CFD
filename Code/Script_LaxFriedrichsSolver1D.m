% Script to test LaxFriedrichsSolver3rOrder1D
% Last update: April 24, 2018
% Solves Solves u_t + f(u)_x = 0    for x in intervalx, t in [0,T]
% with periodic boundary conditions
% Using third order space and time approximations

clear;close all;clc;tic;
% PARAMETERS
f = @(u) u.^2/2;
fp = @(u) u; % derivative of f
intervalx = [-pi pi]; % interval of x
% initial condition
initial_option = 1; % (1) sin x
scheme_option = 4;  % (3) 3rd order
                    % (4) 3rd order with minmod correction
display_option = 2;% (1) error and order of convergence after refinements
                   % (2) plotting an animated solution vs exact solution
                   % (3) plotting refined solutions at the final time step
switch initial_option
    case 1
        u0 = @(x) sin(x);
        U0 = @(x) -cos(x);
        u0p = @(x) cos(x);
end
% discretization parameteres
refinement_size = [40,80,160,320,640];
dx_vector=(intervalx(2)-intervalx(1))./refinement_size;
finalT = 0.5; % the final time
CFL=0.5;
dt_vector=CFL*dx_vector; % time discretization
MAX_ITER = 40; % Max number of iterarions for Newton's method
EPS = 1e-15; % Accuracy of Newton's method

% SOLVER
err_vector=[];
for ind=1:length(dx_vector)
    dx=dx_vector(ind);
    xx = intervalx(1)+dx/2 : dx : intervalx(2)-dx/2;
    dt=dt_vector(ind);
    numFinalT=(floor(finalT/dt))*dt;
    if display_option~=2
        startT=numFinalT-dt;
    else
        startT=0;
    end
    if initial_option == 1
        intervalxP = [intervalx(1)+dx/2 intervalx(2)-dx/2];
        uexacNhalf = BurgersSolver1DPar (u0,u0p,dx,dt,intervalxP,...
            [startT,numFinalT],MAX_ITER,EPS);
        uexacN = BurgersSolver1DPar (u0,u0p,dx,dt,intervalx,...
            [startT,numFinalT],MAX_ITER,EPS);
        % third order numerical integration to approximate cell average
        uexac = (1/6)*(uexacN(:,1:end-1) + uexacN(:,2:end) + 4*uexacNhalf);
    end
    uh = LaxFriedrichsSolver1D(U0,f,fp,dx,dt,intervalx,finalT,scheme_option);
    error=dx*sum(abs(uh(end,:)-uexac(end,:)));
    err_vector=[err_vector;error];%#ok
    if display_option==3
        uh_plot{ind}=uh(end,:);%#ok
    end
end

% DISPLAY
disp('Error and the order of convergence:')
%order = log(err_vector(1:end-1)./err_vector(2:end))/log(2);
order = 0.5*err_vector(1:end-1)./err_vector(2:end);
order = [0;order];
disp([err_vector order]);

switch display_option
    case 2
        disp('Solution for the refined mesh')
        tt = dt : dt : finalT;
        figure
        for nt = 1 : size(uh,1)
            if initial_option == 1
                plot(xx, uh(nt,:),xx,uexac(nt,:));
            end
            title(['time = ' num2str(tt(nt))]);
            pause(0.1);
        end
        legend('Numerical','Exact');
    case 3
        figure;
        for ind=1:length(dx_vector)
            dx=dx_vector(ind);
            xx = intervalx(1)+dx/2 : dx : intervalx(2)-dx/2;
            plot(xx, uh_plot{ind});
            hold on;
            strN = num2str(floor((intervalx(2)-intervalx(1))/dx_vector(ind)));
            legends{ind}=['N = ' strN];%#ok
        end
        title(['Numerical solution using minmod at time T=' num2str(numFinalT)]);
        legend(legends,'Location','southeast');
        xlabel('$x$','interpreter','latex');
        ylabel('Cell averages')
        set(gca,'FontSize',16)    
end
toc;