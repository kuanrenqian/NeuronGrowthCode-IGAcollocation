% Coupled_GreyScott_05202021
% 05/18/2021

%% CleanUp
close all;
clear;
clc;

%% Including Path
addpath('../IGA_collocation_algorithm');

disp('********************************************************************');
disp('2D Reaction Diffusion Gray Scott solver using IGA-Collocation');

%% Variable Initialization
disp('Initializing vairbles...');
% time stepping variables
t = 0;
dt = 0.04;

f = 0.04;
kk = 0.062;
Du  = 1;
Dv  = 0.5;

% k1 = 2;
% k_1 = 10;
% k2 = 20;

% Du = 30;
% Dv = 30;

% diff = 2.5;
% k1b = 0.003;
% k3 = 0.1;
% k4 = 0.5;
% k5 = 1;
% k2 = 20;
% k1 = 0.05;
% alpha1 = 0.6;
% alpha2 = 0.3;
% rhoMax = 2;

% tau = 2;

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
sprs = 1;
[NuNv,N1uNv,NuN1v,~,N2uNv,NuN2v,~,~,size_collpts,~] = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv,sprs);
lap = N2uNv + NuN2v;

%% Iterating variable initialization
disp('Variabel initialization...');
u = ones(size_collpts);
v = zeros(size_collpts);

for i = floor(size_collpts/2)-5:floor(size_collpts/2)+5
    for j = floor(size_collpts/2)-5:floor(size_collpts/2)+5
        u(i,j) = 0;
        v(i,j) = 1;
    end
end

load('phiplot_05202021');
phi = full(phi_plot);
phi= NuNv\reshape(phi,lenu*lenv,1);

% 
% U = rand(size_collpts);
% V = 0.3*rand(size_collpts);


figure(1)
subplot(3,2,1);
contourf(u,'LineStyle','none');
title('Initial - u');
axis square;
colorbar;
subplot(3,2,2);
contourf(v,'LineStyle','none');
title('Initial - v');
axis square;
colorbar;
drawnow;

u = reshape(u,lenu*lenv,1);
v = reshape(v,lenu*lenv,1);
u = NuNv\u;
v = NuNv\v;

u = sparse(u);
v = sparse(v);
lap = sparse(lap);

% Du =sparse(Du);
% Dv =sparse(Dv);

% Ab_atp = zeros(size_collpts);
% Ab_adp = zeros(size_collpts);
% AD = zeros(size_collpts);
% A_adp = zeros(size_collpts);
% D = 0.3*rand(size_collpts);
% A  = rand(size_collpts);
% 
% for i = floor(size_collpts/2)-5:floor(size_collpts/2)+5
%     for j = floor(size_collpts/2)-5:floor(size_collpts/2)+5
%         A(i,j) = 50;
%         D(i,j) = 20;
%     end
% end
% 
% Ab_atp = reshape(Ab_atp,lenu*lenv,1);
% Ab_adp = reshape(zeros(size_collpts),lenu*lenv,1);
% AD = reshape(zeros(size_collpts),lenu*lenv,1);
% A_adp = reshape(zeros(size_collpts),lenu*lenv,1);
% D = reshape(0.3*rand(size_collpts),lenu*lenv,1);
% A  = reshape(rand(size_collpts),lenu*lenv,1);
% 
% Ab_atp_new = Ab_atp;
% Ab_adp_new = Ab_adp;
% AD_new = AD;
% A_new_adp = A_adp;
% D_new = D;
% A_new  = A;

% % ID for boundary location (suppress 4 edges)
% % id = 1 means there is bc
% bcid = zeros([lenu,lenv]);
% for i = 1:lenu
%     bcid(1,i) = 1;
%     bcid(lenu,i) = 1;
%     bcid(i,1) = 1;
%     bcid(i,lenv) = 1;
% end
% bcid = sparse(reshape(bcid,lenu*lenv,1));
% N = zeros(lenu);

for iter = 1:5000
    fprintf('Progress: %.2d\n',iter);
    
%     lapU = (lap*u);
%     lapV = (lap*v);
%     
    %% Grey-Scott model
%     duT = Du.*(lapU)-(NuNv*u).*(NuNv*v).^2.+f.*(1-(NuNv*u));
%     duT = sparse(duT); % 1-NuNv*U is breaking sparse
%     dvT = Dv.*(lapV)+(NuNv*u).*(NuNv*v).^2.-(f+kk).*(NuNv*v);


    %% Grey-Scott
    P = NuNv*phi;
    dPdx = N1uNv*phi;
    dPdy = NuN1v*phi;
    U = NuNv*u;
    dUdx = N1uNv*u;
    dUdy = NuN1v*u;
    lapU = (lap*u);
    V = NuNv*v;
    dVdx = N1uNv*v;
    dVdy = NuN1v*v;
    lapV = (lap*v);    

    %% Uncoupled equations
%     term_Udiff = Du.*(lapU);
%     term_U2 = U.*(V.^2);
%     term_U3 = f.*(1-U);
%     duT = term_Udiff-term_U2+term_U3;
% 
%     term_Vdiff = Dv.*(lapV);
%     term_V2 = U.*(V.^2);
%     term_V3 = (f+kk).*(V);
%     dvT = term_Vdiff+term_V2-term_V3;
    
    %% Coupled equation
    term_Udiff = Du.*(dPdx.*dUdx+dPdy.*dUdy+P.*lapU);
    term_U2 = P.*U.*(V.^2);
    term_U3 = f.*(1-P.*U);
    duT = term_Udiff-term_U2+term_U3;

    term_Vdiff = Dv.*(dPdx.*dVdx+dPdy.*dVdy+P.*lapV);
    term_V2 = U.*(P.*V).^2;
    term_V3 = (f+kk).*(P.*V);
    dvT = term_Vdiff+term_V2-term_V3;
    
    %% simplified activator deactivator model from paper
%     duT = k1.*(NuNv*U).*(1-NuNv*U)-k2.*(NuNv*V).*(t-tau);
%     dvT = D.*(lapV)+k3.*(NuNv*U)-k4.*(NuNv*V);

    %% Actin wave model
%     dAb_atp = (k1.*(NuNv*A)+k1b.*(NuNv*A).*(NuNv*Ab_atp+alpha1.*(NuNv*Ab_adp)+alpha2.*(NuNv*AD))) ...
%         .*((rhoMax-(NuNv*Ab_atp+NuNv*Ab_adp+NuNv*AD))./rhoMax)-k2.*(NuNv*Ab_atp);
%     Ab_atp_new = NuNv\(NuNv*Ab_atp + dAb_atp.*dt);
%     
%     dAb_adp = k2.*(NuNv*Ab_atp)-k3.*(NuNv*Ab_adp).*(NuNv*D);
%     Ab_adp_new = NuNv\(NuNv*Ab_adp + dAb_adp.*dt);
% 
%     dAD = k3.*(NuNv*Ab_adp)-k4.*(NuNv*AD);
%     AD_new = NuNv\(NuNv*AD + dAD.*dt);
%     
%     dA_adp = diff.*(lap*A_adp)+k4.*(NuNv*AD)-k5.*(NuNv*A_adp);
%     A_adp_new = NuNv\(NuNv*A_adp + dA_adp.*dt);
%     
%     dD = diff.*(lap*D)+k4.*(NuNv*AD)-k3.*(NuNv*Ab_adp).*(NuNv*D);
%     D_new = NuNv\(NuNv*D + dD.*dt);    
%     
%     dA = diff.*(lap*A)+k2.*(NuNv*A_adp)-(k1b.*(NuNv*A).*(NuNv*Ab_atp+alpha1.*(NuNv*Ab_adp) ...
%         +alpha2.*(NuNv*AD))+k1.*(NuNv*A)) ...
%         .*(rhoMax-(NuNv*Ab_atp+NuNv*Ab_adp+NuNv*AD)./rhoMax);
%     A_new = NuNv\(NuNv*A + dA.*dt);    


    %% calculating new val
%     % Increments values using Euler's method
%     u_new = NuNv\(NuNv*u + dt*duT);
%     v_new = NuNv\(NuNv*v + dt*dvT);
    u_new = NuNv\(NuNv*u + dt*duT.*P);
    v_new = NuNv\(NuNv*v + dt*dvT.*P);
    
%     [uLHS, uRHS] = StiffMatSetupBCID(NuNv, ((NuNv*U) + dt*duT),bcid,N);
%     [vLHS, vRHS] = StiffMatSetupBCID(NuNv, ((NuNv*V) + dt*dvT),bcid,N);
%     U_new = uLHS\uRHS;
%     V_new = vLHS\vRHS;
   
%% Plotting figures
    if(mod(iter,25) == 0 )
        u_plot = reshape(NuNv*u_new,size_collpts,size_collpts);
        subplot(3,2,3);
        imagesc(u_plot);
        title(sprintf('u at time = %.2d',t));
        axis square;
        colorbar;

        v_plot = reshape(NuNv*v_new,size_collpts,size_collpts);
        subplot(3,2,4);
        imagesc(v_plot);
        title(sprintf('v at time = %.2d',t));
        axis square;
        colorbar;

        subplot(3,2,5);
        imagesc(reshape(NuNv*phi,lenu,lenv));
        title(sprintf('phi at time = %.2d',t));
        axis square;
        colorbar;

    %     A_plot = reshape(NuNv*A_new,size_collpts,size_collpts);
    %     subplot(1,2,1);
    %     imagesc(A_plot);
    %     title(sprintf('A at time = %.2d',t));
    %     axis square;
    %     colorbar;
    %     
    %     D_plot = reshape(NuNv*D_new,size_collpts,size_collpts);
    %     subplot(1,2,2);
    %     imagesc(D_plot);
    %     title(sprintf('D at time = %.2d',t));
    %     axis square;
    %     colorbar;

        % plot current iteration
        drawnow;
    end

    %% iteration update
    % update variables in this iteration    
    u = u_new;
    v = v_new;
% 
%     Ab_atp = Ab_atp_new;
%     Ab_adp = Ab_adp_new;
%     AD = AD_new;
%     A_adp = A_new_adp;
%     D = D_new;
%     A = A_new;

    t = t+dt;
    
end
disp('All simulations complete!\n');
disp('********************************************************************');
