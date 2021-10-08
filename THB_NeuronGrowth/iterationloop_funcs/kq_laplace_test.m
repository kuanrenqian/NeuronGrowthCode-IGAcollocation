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

% getting active coordinates
ac_len = length(ac);
ac_x = zeros(ac_len,1);
ac_y = zeros(ac_len,1);
for i = 1:ac_len
    ac_level = ac(i,2);
    ac_indx = ac(i,1);
    ac_x(i) = cell2mat(Em{ac_level}(ac_indx,8));
    ac_y(i) = cell2mat(Em{ac_level}(ac_indx,9));
end

% setting bcid (dirichlet boundary condition id)
bc_layers = 2;
[bcid,bcid_cp] = kqMakeBCID(coll_sz,bc_layers,coll_X,coll_Y,THBfinal);

T_initial_cp = T_analytical_cp;
T_initial_cp(bcid_cp==0) = 0;

% solving linear sustem
[lhs, rhs] = kqLaplaceStiffMatSetupBCID(cm.lap,T_initial_cp,bcid_cp);

T_sol = lhs\rhs;

% interpolate from locally refined points to structured grid for plotting
x_sol = THBfinal(bcid_cp==0,1);
y_sol = THBfinal(bcid_cp==0,2);
T_sol_plot  = griddata(x_sol,y_sol,T_sol,coll_X,coll_Y);
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
