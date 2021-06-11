% Isogeometric-Collocation method (IGA-Collocation) for Gray Scott model on 2D plate
% Kuanren Qian
% 03/11/2021

%% CleanUp
close all;
clear;
clc;

%% Including Path
addpath('../IGA_collocation_algorithm');
addpath('../github_repo');

disp('********************************************************************');
disp('2D Reaction Diffusion Gray Scott solver using IGA-Collocation');

%% Variable Initialization
disp('Initializing vairbles...');
% time stepping variables
scale = 30;
t = 0;
dt = 1/scale;

end_t = 2000;

f  = 0.04;
kk  = 0.0636;
Du  = 1;
Dv  = 0.5;
% 
% f = f/scale;
% kk = kk/scale;
% Du = Du/scale;
% Dv = Dv/scale;

% B-spline curve order (U,V direction)
p = 3;
q = 3;

% Control mesh size (knot vector)
Nx = 60;
Ny = 60;
dx = 1;
dy =1;
knotvectorU = [zeros([1,p]),0:Nx,ones([1,p])*Nx].';
knotvectorV = [zeros([1,p]),0:Ny,ones([1,p])*Ny].';
lenu = length(knotvectorU)-2*(p-1);
lenv = length(knotvectorV)-2*(p-1);

%% Constructing coef matrix (concise function version)
order_deriv = 2;
sprs = 0;
[NuNv,~,~,~,N2uNv,NuN2v,~,~,size_collpts,~] = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv,sprs);
lap = N2uNv + NuN2v;

%% Iterating variable initialization
disp('Variabel initialization...');
U = ones(size_collpts);
V = zeros(size_collpts);

for i = floor(size_collpts/2)-5:floor(size_collpts/2)+5
    for j = floor(size_collpts/2)-5:floor(size_collpts/2)+5
        U(i,j) = 0;
        V(i,j) = 1;
    end
end

figure(1)
subplot(1,2,1);
contourf(U,'LineStyle','none');
title('Initial - u');
axis square;
colorbar;
subplot(1,2,2);
contourf(V,'LineStyle','none');
title('Initial - v');
axis square;
colorbar;
drawnow;

U = reshape(U,lenu*lenv,1);
V = reshape(V,lenu*lenv,1);
% 
% U = sparse(U);
% V = sparse(V);
% lap = sparse(lap);
% NuNv = sparse(NuNv);

while t<end_t    
    fprintf('Progress: %.2d/%.2d\n',t,end_t);
    
    lapU = (lap*U)/scale;
    lapV = (lap*V)/scale;

%     lapU = reshape(laplacian(reshape((NuNv*U),lenu,lenv)),lenu*lenv,1);
%     lapV = reshape(laplacian(reshape((NuNv*V),lenu,lenv)),lenu*lenv,1);
    
    duT = (Du.*(lapU)-(NuNv*U).*(NuNv*V).^2.+f.*(1-(NuNv*U)));
    dvT = (Dv.*(lapV)+(NuNv*U).*(NuNv*V).^2.-(f+kk).*(NuNv*V));

    
    duT = duT*50;
    dvT = dvT*50;

    % Increments values using Euler's method
    U_new = NuNv\((NuNv*U) + dt*duT);
    V_new = NuNv\((NuNv*V) + dt*dvT);
   
%% Plotting figures
    u_plot = reshape(NuNv*U_new,size_collpts,size_collpts);
    subplot(1,2,1);
    imagesc(u_plot);
    title(sprintf('u at time = %.2d',t));
    axis square;
    colorbar;
    
    v_plot = reshape(NuNv*V_new,size_collpts,size_collpts);
    subplot(1,2,2);
    imagesc(v_plot);
    title(sprintf('v at time = %.2d',t));
    axis square;
    colorbar;

    % plot current iteration
    drawnow;
    
    %% iteration update
    % update variables in this iteration    
    U = U_new;
    V = V_new;
    t = t+dt;
    
end
disp('All simulations complete!\n');
disp('********************************************************************');
