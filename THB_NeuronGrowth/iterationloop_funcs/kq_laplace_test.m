%% Solving laplace on THB splines
% calculating THB NuNv,N1uNv... and store in cm
cm = kqTHBders(Pixel,Jm,Pm,Dm);

% Get locally refined points from THB
THBfinal = zeros(bf_ct,2);
for i = 1:bf_ct
    bbc = bf(i,1:2);
    bf_lev = bf(i,3);
    bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
    Pi = Pm{bf_lev,1};
    THBfinal(i,1) = Pi(bi,1);
    THBfinal(i,2) = Pi(bi,2);
end

% Calculating analytical solution on structured collocation grid, then
% interp to locally refined mesh 
CEb = Pm{1,1};
[M,~] = size(CEb);
coll_sz = sqrt(M);

[M,~] = size(THBfinal);
T_analytical_cp = zeros(M,1);
for k = 1:M
    x = THBfinal(k,1);
    y = THBfinal(k,2);
    T_analytical_cp(k) = x*y+x+y+1;
%     T_analytical(k) = 2*x^2+2*y^2+x*y-x-y+1;
end

% create collocation mesh grid
coll_X_array = CEb(1:coll_sz,2).';
coll_Y_array = CEb(1:coll_sz,2).';
[coll_X,coll_Y] = meshgrid(coll_X_array,coll_Y_array);

% setting bcid (dirichlet boundary condition id)
bc_layers = 1;
[bcid,bcid_cp] = kqMakeBCID(coll_sz,bc_layers,coll_X,coll_Y,THBfinal);

T_initial_cp = T_analytical_cp;
T_initial_cp(bcid_cp==0) = 0;

% solving linear sustem
RHS = -cm.lap*(T_initial_cp);
[LHS, c] = kqStiffMatSetupBCID(cm.lap,RHS,bcid_cp);
T_sol = LHS\(RHS);

% interpolate from locally refined points to structured grid for plotting
T_sol_plot  = griddata(THBfinal(:,1),THBfinal(:,2),T_sol,coll_X,coll_Y);
T_ana_plot  = griddata(THBfinal(:,1),THBfinal(:,2),T_analytical_cp,coll_X,coll_Y);

%% Plot all figures
figure;
colormap parula;
subplot(2,2,1);
imagesc(T_ana_plot(bc_layers+3:end-bc_layers-2,bc_layers+3:end-bc_layers-2));
title('Analytical Laplace Solution');
colorbar;
axis square;

subplot(2,2,2);
% if(maxlev>1)
%     imagesc([0,1],[0,1],reshape(grad_log,coll_sz,coll_sz)); 
%     hold on;
% end
displayAdaptiveGrid(ac,Coeff,Em,knotvectorU,knotvectorV,Jm,Pm,parameters,dx*nx,dy*ny);
title(sprintf('Mesh with %d refinements',maxlev-1));
axis square;

subplot(2,2,3);
T_sol_plot = reshape(T_sol_plot,coll_sz,coll_sz);
imagesc(T_sol_plot(bc_layers+3:end-bc_layers-2,bc_layers+3:end-bc_layers-2));
title('Calculated solution');
colorbar;
axis square;

subplot(2,2,4);
% extracting points (avoiding boundary points)
T_ana = T_ana_plot(bcid==0);
T_solu = T_sol_plot(bcid==0);
%calculate absolute error on extracted points
T_abs_err = reshape(T_ana-T_solu,coll_sz-(bc_layers+1)*2,coll_sz-(bc_layers+1)*2);
imagesc(T_abs_err(2:end-1,2:end-1));
title('Absolute Error');
colorbar;
axis square;
drawnow;   

%%
% figure;
% scatter3(THBfinal(:,1),THBfinal(:,2),T_sol,140,T_sol,'fill','Marker', 's');
% view(90,90)
% colorbar;

%%

% [M,~] = size(THBfinal);
% T_analytical_coll = zeros(M,1);
% for k = 1:M
%     x = THBfinal(k,1);
%     y = THBfinal(k,2);
%     T_analytical_coll(k) = x*y+x+y+1;
% %     T_analytical(k) = 2*x^2+2*y^2+x*y-x-y+1;
% end
% 
% % cp -> control points, coll -> collocation points
% T_analytical_cp = cm.NuNv\T_analytical_coll;
% % T_analytical_coll_1 = cm.NuNv*T_analytical_cp;
% % T_ana_plot  = griddata(THBfinal(:,1),THBfinal(:,2),T_analytical_coll_1,coll_X,coll_Y);
% 
% % figure
% % plotp = reshape(T_analytical_coll,22,22)-reshape(T_analytical_coll_1,22,22);
% % imagesc(plotp(5:17,5:17));
% 
% % create collocation mesh grid
% cp_X_array = linspace(0,1,Nx).';
% cp_Y_array = linspace(0,1,Nx).';
% [cp_X,cp_Y] = meshgrid(cp_X_array,cp_Y_array);
% 
% % setting bcid (dirichlet boundary condition id)
% bc_layers = 1;
% 
% [bcid,bcid_cp] = kqMakeBCID(Nx,bc_layers,cp_X,cp_Y,THBfinal);
% 
% T_initial_cp = T_analytical_cp;
% T_initial_cp(bcid==0) = 0;
% 
% % solving linear sustem
% RHS = -cm.lap*(T_initial_cp);
% [LHS, c] = kqStiffMatSetupBCID(cm.lap,RHS,bcid);
% T_sol = LHS\(RHS);
% 
% % % interpolate from locally refined points to structured grid for plotting
% % T_S  = scatteredInterpolant(THBfinal(:,1),THBfinal(:,2),cm.NuNv*T_sol);
% % T_sol_plot = reshape(T_S.Values,coll_sz,coll_sz);
% % T_A  = scatteredInterpolant(THBfinal(:,1),THBfinal(:,2),T_analytical_coll);
% % T_ana_plot = reshape(T_A.Values,coll_sz,coll_sz);
% 
% % interpolate from locally refined points to structured grid for plotting
% T_sol_plot  = griddata(THBfinal(:,1),THBfinal(:,2),cm.NuNv*T_sol,coll_X,coll_Y);
% T_ana_plot  = griddata(THBfinal(:,1),THBfinal(:,2),T_analytical_coll,coll_X,coll_Y);
% 
% figure;
% subplot(2,2,1);
% imagesc(T_ana_plot);
% title('Analytical solution');
% colorbar;
% 
% subplot(2,2,2);
% displayAdaptiveGrid(ac,Coeff,Em,knotvectorU,knotvectorV,Jm,Pm,parameters,dx*nx,dy*ny);
% title(sprintf('Mesh with %d refinements',maxlev-1));
% axis square;
% 
% subplot(2,2,3);
% imagesc(T_sol_plot);
% title('Calculated solution');
% colorbar;
% 
% subplot(2,2,4);
% imagesc(T_sol_plot-T_ana_plot);
% title('Absolute error');
% colorbar;
% 
% %%
% % phi1 = reshape(phi,484,1);
% % pp1 = cm.NuNv\phi1;
% % pp2 = cm.NuNv*pp1;
% % 
% % figure;
% % subplot(1,3,1);
% % plot1 = reshape(phi,22,22);
% % imagesc(plot1);
% % subplot(1,3,2);
% % plot2 = reshape(pp1,20,20);
% % imagesc(plot2);
% % subplot(1,3,3);
% % plot3 = reshape(pp2,22,22);
% % imagesc(plot3);
