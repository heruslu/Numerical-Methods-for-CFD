% Script to test Burgers Solver
% Last update: March 14, 2018
close all; clear; clc; tic;
% initial condition
u0 = @(x) sin(x);
u0p = @(x) cos(x);
figure_option = 2;
% 1 for plot of solutions at some given time steps
% 2 for plotting an animated solution
% discretization parameteres
intervalx = [-7 7];
dx=0.1; % space discretization
xx = intervalx(1) : dx : intervalx(2);
Tstar = 1; % the time where shock developes
intervalt = [0 Tstar];
dt=0.01; % time discretization
MAX_ITER = 20; % Max number of iterarions for Newton's method which occurs
EPS = 1e-14; % Accuracy of Newton's method
% in the solver
u = BurgersSolver1DPar (u0,u0p,dx,dt,intervalx,intervalt,MAX_ITER,EPS);
switch figure_option
    case 1
        time_vec = [Tstar/3, Tstar/2, 2*Tstar/3, 4*Tstar/5];
        time_ind = floor(time_vec/dt);
        plot(xx, u(time_ind(1),:),xx, u(time_ind(2),:),xx, ...
            u(time_ind(3),:),xx, u(time_ind(4),:));
        title('Burgers Solution','interpreter','latex')
        str1 = '$u(x,1/3)$'; str2 = '$u(x,1/2)$'; str3 = '$u(x,2/3)$';
        str4 = '$u(x,4/5)$';
        h = legend(str1,str2,str3,str4,'Location','northeast');
        set(h,'interpreter','latex')
        set(gca,'FontSize',16);
        xlim(intervalx);
    case 2
        tt = dt : dt : Tstar;
        Nt = size(tt,2);
        figure
        for nt = 1 : Nt
            plot(xx, u(nt,:));
            title(['time = ' num2str(tt(nt))]);
            pause(0.01);
        end
end
toc;