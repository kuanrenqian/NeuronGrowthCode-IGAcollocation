% IGA-collocation Implementation with laplace and poisson problems
% Kuanren Qian
% 04/20/2021

%% CleanUp
close all;
clear;
clc;

%% Including Path
addpath('../IGA_collocation_algorithm');

disp('**************************************************************************');
disp('IGA-collocation Implementation with laplace and poisson problems');
disp('Kuanren Qian - 04/20/2021');

%% Variable Initialization
disp('Initializing vairbles...');
% time stepping variables (60s transient simulation)
t = 0;
dt = 0.01;
end_t = 40;

% heat transfer variables (Copper)
k_heat = 385.0; %W/mK
rho = 8960; %kg/m^3
Cp = 0.39; %kJ/kgK
alpha = k_heat/(rho*Cp);

% B-spline curve order (U,V direction)
p = 3;
q = 3;

% for L2 convergence study
grid_size_X = [4,8,16,32,64];
grid_size_Y = [4,8,16,32,64];
L2norm = zeros([5,1]);

grid_size_X = 10;
grid_size_Y = 10;

syms X Y;
% pool = parpool;

% looping through grid sizes for convergence study
for n = 1:length(grid_size_X)
    tic
    % Control mesh size (knot vector)
    grid_size_x = grid_size_X(n);
    grid_size_y = grid_size_Y(n);
    knotvectorU = [0,0,0,linspace(0,1,grid_size_x),1,1,1].';
    knotvectorV = [0,0,0,linspace(0,1,grid_size_y),1,1,1].';
    
    % setting lenu lenv this way for easier access to ghost nodes later on
    % (1 layer of ghost nodes all around)
    lenu = length(knotvectorU)-2*(p-1);
    lenv = length(knotvectorV)-2*(q-1);

    %% Constructing coef matrix (concise function version)
    order_deriv = 2;
    sprs = 1;
    [NuNv,N1uNv,NuN1v,N1uN1v,N2uNv,NuN2v,N2uN2v,coll_p,size_collpts,Control_points] = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv, sprs);
    lap = N2uNv + NuN2v;

    %% Steady state computation
    fprintf('Steady State - DOF: %d\n', grid_size_x*grid_size_y);
    T_analytical = zeros([lenu,lenv]);
    parfor i = 1:lenu
        for j = 1:lenv
            k = (i-1)*lenu+j;
            x = coll_p(k,1);
            y = coll_p(k,2);
            %% Laplace - Xiaodong analytical solution
%             theta = atan(y/x);
%             r = x/cos(theta);
%             T_analytical(i,j) = r^(2/3)*sin(2*theta/3-pi/3);
            %% Laplace - simple analytical solution
%             T_analytical(i,j) = x*y+x+y+1;
            %% Laplace - sinxcoshy analytical solution
%             T_analytical(i,j) = sin(x)*cosh(y);
            %% Laplace - sin2xsinh2y analytical solution
%             T_analytical(i,j) = sin(2*x)*sinh(2*y);
            %% Poisson - simple analytical solution
%             T_analytical(i,j) = 2*x^2+2*y^2+x*y-x-y+1;
            %% Poisson - convergence paper analytical solution
            % https://www.researchgate.net/publication/301873952_Isogeometric_Least-Squares_Collocation_Method_with_Consistency_and_Convergence_Analysis
            T_analytical(i,j) = (x^2+y^2-1)*(x^2+y^2-16)*sin(x)*sin(y);
            f(i,j) = (3*x^4-67*x^2-67*y^2+3*y^4+6*x^2*y^2+116)*sin(x)*sin(y)...
                +(68*x-8*x^3-8*x*y^2)*cos(x)*sin(y)...
                +(68*y-8*y^3-8*y*x^2)*cos(y)*sin(x);
           %% Poisson - NIST paper analytical solution (NIST-1)
%             % https://core.ac.uk/download/pdf/82459498.pdf
%             T_analytical(i,j) = 2^40*x^10*(1-x)^10*y^10*(1-y)^10;
%             F_analytical = @(X,Y) 2^40*X^10*(1-X)^10*Y^10*(1-Y)^10;
%             D2Fdx2 = diff(F_analytical,X,2);
%             D2Fdy2 = diff(F_analytical,Y,2);
%             d2fdx2 = subs(D2Fdx2,X,x);
%             d2fdx2 = subs(d2fdx2,Y,y);
%             d2fdy2 = subs(D2Fdy2,X,x);
%             d2fdy2 = subs(d2fdy2,Y,y);
%             f(i,j) = -double(d2fdx2 + d2fdy2);

            %% Poisson - NIST paper analytical solution (NIST-4)
%             % https://core.ac.uk/download/pdf/82459498.pdf
%             x_loc = 0.5;
%             y_loc = 0.5;
%             alph = 1000;
%             T_analytical(i,j) = exp(-alph*((x-x_loc)^2+(y-y_loc)^2));
%             F_analytical = @(X,Y) exp(-alph*((X-x_loc)^2+(Y-y_loc)^2));
%             D2Fdx2 = diff(F_analytical,X,2);
%             D2Fdy2 = diff(F_analytical,Y,2);
%             d2fdx2 = subs(D2Fdx2,X,x);
%             d2fdx2 = subs(d2fdx2,Y,y);
%             d2fdy2 = subs(D2Fdy2,X,x);
%             d2fdy2 = subs(d2fdy2,Y,y);
%             f(i,j) = -double(d2fdx2 + d2fdy2);
            
        end
    end
    
    %% BC assignment (Dirichlet BC all around)
    T_bc = NuNv\reshape(T_analytical,lenu*lenv,1);
    T_initial = reshape(T_bc,lenu,lenv);
    for i = 2:length(T_initial)-1
        for j = 2:length(T_initial)-1
            T_initial(i,j) = 0;
        end
    end
    T_initial = reshape(T_initial,(lenu)^2,1);

    % ID for non-zero boundary location
    % id = 1 means there is bc
    bcid = zeros([lenu,lenv]);
    for i = 1:lenu
        bcid(1,i) = 1;
        bcid(lenu,i) = 1;
        bcid(i,1) = 1;
        bcid(i,lenv) = 1;
    end
    bcid = reshape(bcid,lenu*lenv,1);
    
    %% Solve Laplace
%     lap_ss = N2uNv+NuN2v;
%     b = -lap*(T_initial);
%     [coll_Lhs, coll_Rhs] = StiffMatSetupBCID(coll_Lhs, coll_Rhs,bcid,T_initial);

%     T_new = laplacian_ss\(b);

    %% Solve Simple Poisson
%     coll_lhs = N2uNv+NuN2v;
%     source_term = 8*(NuNv*ones(size(T_initial)));
%     coll_rhs = -lap*(T_initial)+source_term;
%     [coll_Lhs, coll_Rhs] = StiffMatSetupBCID(coll_Lhs, coll_Rhs,bcid,T_initial);

%     T_new = coll_lhs\(coll_rhs);

    %% Solve Poisson from convergence paper
    coll_Lhs = -lap+NuNv;
    coll_Rhs = reshape(f,lenu*lenv,1)-coll_Lhs*T_initial;
    [coll_Lhs, coll_Rhs] = StiffMatSetupBCID(coll_Lhs, coll_Rhs,bcid,T_initial);
    
    T_new = coll_Lhs\(coll_Rhs);
    
    %% Solve Poisson from NIST paper (NIST-1&4)
%     coll_Lhs = -lap;
%     coll_Rhs = reshape(f,lenu*lenv,1)-coll_Lhs*T_initial;
%     [coll_Lhs, coll_Rhs] = StiffMatSetupBCID(coll_Lhs, coll_Rhs,bcid,T_initial);
%     
%     T_new = coll_Lhs\(coll_Rhs);
    
    %% calculating and recording L2 norm error convergence
    T_p =  reshape(NuNv*T_new,lenu,lenv);
    L2error = kqL2normError(T_analytical,T_p);
    L2norm(n) = L2error;
    fprintf('L2 norm error: %f\n', L2error);
    toc
end

%% Plotting Results and error

% figure(1);
% plot(Control_points(:,1),Control_points(:,2),'xb',coll_p(:,1),coll_p(:,2),'or');
% legend('Control points - knot vector','Collocation Points');
% title('Control Mesh and Collocation Points');
% disp('Done!');

% plotting analytical solution
set(gcf,'position',[100,100,1500,300])
subplot(1,4,1);
imagesc(T_analytical);
title('Analytical Solution');
xlabel('Grid size (coord 0~2)');
ylabel('Grid size (coord 0~2)');
colormap (gca,'parula');
colorbar;

% plotting calculated IGA-collocation results
subplot(1,4,2);
T_plot = reshape(NuNv*T_new,lenu,lenv);
imagesc(T_plot(2:end-1,2:end-1));
title('IGA-collocation results');
xlabel('Grid size (coord 0~2)');
ylabel('Grid size (coord 0~2)');
colormap (gca,'parula');
colorbar;

% plotting error distribution
subplot(1,4,3);
error_plot = (T_analytical(2:end-1,2:end-1)-T_p(2:end-1,2:end-1));
imagesc(error_plot);
title('Absolute Error');
xlabel('Grid size (coord 0~2)');
ylabel('Grid size (coord 0~2)');
colormap (gca,'jet');
colorbar;

dof = grid_size_X(1:n).*grid_size_Y(1:n);
subplot(1,4,4);
loglog(dof,L2norm);
title('L2 norm error');
xlabel('DOF');
ylabel('L2 norm error');
grid on;
% legend(problem);

% logx = log(dof);
% logy = log(L2norm);
% Const = polyfit(logx, logy, 1);
% hold on
% plot(dof, exp(polyval(Const, logx)),'c*');

% poolobj = gcp('nocreate');
% delete(poolobj);

disp('Simulation Completed!');