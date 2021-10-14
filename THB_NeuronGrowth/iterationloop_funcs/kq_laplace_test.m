%% Solving laplace on THB splines
% calculating THB NuNv,N1uNv... and store in cm
cm = kqTHBders(Pixel,Jm,Pm,Dm);

% Calculating analytical solution on locally refined mesh 
M = length(THBfinal);
T_analytical_cp = zeros(M,1);
for k = 1:M
    x = THBfinal(k,1);
    y = THBfinal(k,2);
    T_analytical_cp(k) = x*y+x+y+1;
%     T_analytical(k) = 2*x^2+2*y^2+x*y-x-y+1;
end

%%
% setting bcid (dirichlet boundary condition id)
bc_layers = 2;
[bcid,bcid_cp] = kqMakeBCID(Nx,bc_layers,X,Y,THBfinal);

T_initial_cp = T_analytical_cp;
T_initial_cp(bcid_cp==0) = 0;

% solving linear sustem
[lhs, rhs] = kqLaplaceStiffMatSetupBCID(cm.lap,T_initial_cp,bcid_cp);

T_sol = lhs\rhs;

% interpolate from locally refined points to structured grid for plotting
x_sol = THBfinal(bcid_cp==0,1);
y_sol = THBfinal(bcid_cp==0,2);
T_sol_plot  = griddata(x_sol,y_sol,T_sol,X,Y);
T_ana_plot  = griddata(THBfinal(:,1),THBfinal(:,2),T_analytical_cp,X,Y);

T_ana_plot(bcid==1)=[];
T_sol_plot(bcid==1)=[];
sz = sqrt(length(T_ana_plot));
T_ana_plot = reshape(T_ana_plot,sz,sz);
T_sol_plot = reshape(T_sol_plot,sz,sz);

%% Plot all figures
figure;
colormap parula;
subplot(2,2,1);
imagesc(T_ana_plot);
title('Analytical Laplace Solution');
colorbar;
axis square;

subplot(2,2,2);
displayAdaptiveGrid(ac,Coeff,Em,knotvectorU,knotvectorV,Jm,Pm,parameters);
title(sprintf('Mesh with %d refinements',maxlev-1));
axis square;

subplot(2,2,3);
imagesc(T_sol_plot);
title('Calculated solution');
colorbar;
axis square;

subplot(2,2,4);
% calculate absolute error on extracted points
T_abs_err = T_ana_plot - T_sol_plot;
imagesc(T_abs_err);
title('Absolute Error');
colorbar;
axis square;
drawnow;