% IGA-collocation Implementation for Transient 2D heat transfer
% Kuanren Qian
% 04/20/2021

%% CleanUp
close all;
clear;
clc;

format LONG;

%% Including Path
addpath('../IGA_collocation_algorithm');

disp('********************************************************************');
disp('IGA-collocation Implementation for Transient 2D heat transfer');
disp('Kuanren Qian - 04/22/2021');

%% Variable Initialization
disp('Initializing vairbles...');
% time stepping variables (60s transient simulation)
t = 0;
dt = 0.001;
end_t = 40;

% heat transfer variables (Copper)
k_heat = 385.0; %W/mK
rho = 8960; %kg/m^3
Cp = 0.39; %kJ/kgK
alpha = k_heat/(rho*Cp);

% B-spline curve order (U,V direction)
p = 3;
q = 3;

% looping through grid sizes for convergence study
% Control mesh size (knot vector)
grid_size_x = 60;
grid_size_y = 60;
knotvectorU = [zeros([1,p]),0:grid_size_x,ones([1,p])*grid_size_x].';
knotvectorV = [zeros([1,q]),0:grid_size_y,ones([1,q])*grid_size_y].';

% setting lenu lenv this way for easier access to ghost nodes later on
% (1 layer of ghost nodes all around)
lenu = length(knotvectorU)-2*(p-1);
lenv = length(knotvectorV)-2*(q-1);

%% Constructing coef matrix (concise function version)
order_deriv = 2;
sprs = 0;
[NuNv,N1uNv,NuN1v,N1uN1v,N2uNv,NuN2v,N2uN2v,coll_p,size_collpts,Control_points] = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv,sprs);
lap = N2uNv + NuN2v;

set(gcf,'position',[100,100,800,300]);
subplot(1,5,1);
plot(Control_points(:,1),Control_points(:,2),'xb',coll_p(:,1),coll_p(:,2),'or');
legend('Control points - knot vector','Collocation Points');
title('Control Mesh and Collocation Points');
disp('Done!');

%% Transient
disp('Initializing temperature field...');
T_all = zeros(lenu);
for i = 1:floor(length(T_all)/2)
    T_all(i,1:2) = [1,1];
end
T_all = reshape(T_all,(lenu)^2,1);
T_initial = T_all;

% ID for non-zero boundary location
% id = 1 means there is bc
bcid = zeros([lenu,lenv]);
for i = 1:lenu
    bcid(i,1) = 1;
    bcid(i,2) = 1;
end
bcid = reshape(bcid,lenu*lenv,1);

% plotting initial temperature field
subplot(1,2,1);
T_initial_plot = reshape(T_initial,lenu,lenv);
contourf(T_initial_plot,20,'LineStyle','none');
axis square;
title('Initial Temperature');
colorbar;
disp('Done!');

disp('Transient iterations...');
while t<=end_t
    disp(t);
    % Calculate Temperature in the next time step
    coll_Lhs = NuNv-dt*alpha*(lap);
    coll_Rhs = NuNv*T_all;
    
    [coll_Lhs, coll_Rhs] = StiffMatSetupBCID(coll_Lhs, coll_Rhs,bcid,NuNv*T_initial);

    T_new = coll_Lhs\coll_Rhs;

    % time increment
    t = t+dt;

    T_plot = reshape(NuNv*T_new,lenu,lenv);
    subplot(1,2,2);
    imagesc(T_plot);
    contourf(T_plot,20,'LineStyle','none');
    axis square;
    title(sprintf('Temperature (K) at t=%.2fs',t));
    colorbar;
    drawnow;

    T_all = T_new;
end
disp('Done!');

