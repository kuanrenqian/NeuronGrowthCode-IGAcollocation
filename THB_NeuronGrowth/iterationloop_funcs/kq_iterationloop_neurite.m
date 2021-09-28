%% Solving laplace on THB splines
% calculating THB NuNv,N1uNv... and store in cm
cm = kqTHBders(Pixel,Jm,Pm,Dm);

CEb = Pm{1,1};
[M,~] = size(CEb);
coll_sz = sqrt(M);

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
% create collocation mesh grid
coll_X_array = CEb(1:coll_sz,2).';
coll_Y_array = CEb(1:coll_sz,2).';
[coll_X,coll_Y] = meshgrid(coll_X_array,coll_Y_array);
cp_X_array = linspace(0,1,Nx).';
cp_Y_array = linspace(0,1,Nx).';
[cp_X,cp_Y] = meshgrid(cp_X_array,cp_Y_array);

% setting bcid (dirichlet boundary condition id)
bc_layers = 1;
[bcid,bcid_cp] = kqMakeBCID(Nx,bc_layers,cp_X,cp_Y,THBfinal);

% % % % initializing phi, theta, and temperature   
% len_cp = length(THBfinal);
% phi_cp = zeros(len_cp,1);
% conct_cp = zeros(len_cp,1);
% for i = 1:len_cp
%     x = THBfinal(i,1);
%     y = THBfinal(i,2);
%     if (sqrt((x-0.5)^2+(y-0.5)^2) < seed_radius*dx)
%         r = sqrt((x-0.5)^2+(y-0.5)^2);
%         phi_cp(i) = 1;
%         conct_cp(i) = (0.5+0.5*tanh(seed_radius-r)/2);
%     end
% end
% theta_cp  = rand([len_cp,1]);
% tempr_cp = zeros(len_cp,1);
% theta_ori_cp = zeros(len_cp,1);

phi_cp = cm.NuNv\phi;
theta_cp = cm.NuNv\rand(Nx*Ny,1);
tempr_cp = cm.NuNv\zeros(Nx*Ny,1);

% Visualize initialization of phi on locally refined collocation points
figure
subplot(2,2,1);
% scatter3(THBfinal(:,1),THBfinal(:,2),phi_cp,30,phi_cp,'fill');
imagesc(reshape(cm.NuNv*phi_cp,Nx,Ny));
title('Initial Phi');
subplot(2,2,2);
% scatter3(THBfinal(:,1),THBfinal(:,2),theta_cp,30,theta_cp,'fill');
imagesc(reshape(cm.NuNv*theta_cp,Nx,Ny));
title('Initial Theta');
subplot(2,2,3);
% scatter3(THBfinal(:,1),THBfinal(:,2),tempr_cp,30,tempr_cp,'fill');
imagesc(reshape(cm.NuNv*tempr_cp,Nx,Ny));
title('Initial T');
% subplot(2,2,4);
% scatter3(THBfinal(:,1),THBfinal(:,2),conct_cp,30,conct_cp,'fill');
% title('Initial Tubulin');

phi_initial = zeros(size(phi_cp));
% conct_initial  = zeros(size(conct_cp));
theta_initial  = theta_cp;
tempr_initial  = zeros(size(tempr_cp));

disp('Phi,conc,theta,tempr,bcid - initialization done!');
disp('********************************************************************');

kq_growth_stage_var_initialization

figure;
tic
% transient iterations
for iter=1:1:end_iter

%     if(mod(iter,10) == 0)
        fprintf('Progress: %.2d/%.2d\n',iter,end_iter);
        toc
%     end

    % calculating a and a*a' (aap) in the equation using theta and phi
    [a_cp, ~, aap_cp,~,~] = kqGetEpsilonAndAap(epsilonb,delta,phi_cp,theta_cp,cm.NuNv,cm.NuN1v,cm.N1uNv);

    tempr = cm.NuNv*tempr_cp;
    E = (alph./pix).*atan(gamma.*(1-tempr));

%     if(iter<=iter_stage2_begin)
%         E = (alph./pix).*atan(gamma.*(1-tempr));
%     else
%         conct = cm.NuNv*conct_cp;
%         r = 1e6;
%         g = 0.1;
%         delta_L = r*conct - g;
%         term_change = (regular_Heiviside_fun(delta_L));
%         term_change(theta_ori==1)=1;
%         E = (alph./pix).*atan(term_change.*gamma.*(1-tempr));
%         E(abs(theta_ori)==0) = 0;
%         
%         if(mod(iter,50) == 0)
%             subplot(2,2,2);
%             phi_plot  = griddata(THBfinal(:,1),THBfinal(:,2),phi_cp,cp_X,cp_Y);
% 
%             imagesc(E+phi_plot);
%             title(sprintf('E overlay with phi'));
%             axis square;
%             colorbar;
%         end
%     end

    %% Phi (Implicit Nonlinear NR method)
    % NR method initial guess (guess current phi)
    phiK_cp = phi_cp;
    NNp = cm.NuNv*phi_cp;

    % residual for NR method (make sure it is larger than tolerance at the beginning)
    R = 2*tol;

    % magnitude of theta gradient (offset by 1e-12 to avoid division by zero)
    mag_grad_theta = sparse(sqrt((cm.N1uNv*theta_cp).^2+(cm.NuN1v*theta_cp).^2)+1e-12);
    C1 = E-0.5+6.*s_coeff.*mag_grad_theta;

    % NR method calculation
    ind_check = 0;
    NNa = cm.NuNv*a_cp;
    N1Na = cm.N1uNv*a_cp;
    NN1a = cm.NuN1v*a_cp;
    NNaap = cm.NuNv*aap_cp;
    N1Naap = cm.N1uNv*aap_cp;
    NN1aap = cm.NuN1v*aap_cp;

    out = cell(4,1);
    dt_t = 0;
    while max(abs(R)) >= tol
        NNpk = cm.NuNv*phiK_cp;
        N1Npk = cm.N1uNv*phiK_cp;
        NN1pk = cm.NuN1v*phiK_cp;
        N1N1pk = cm.N1uN1v*phiK_cp;
        LAPpk = cm.lap*phiK_cp;

        % term a2
        out{1}  = 2*NNa.*N1Na.*N1Npk+NNa.^2.*LAPpk ...
            +2*NNa.*NN1a.*NN1pk;
        % termadx
        out{2} = N1Naap.*NN1pk+NNaap.*N1N1pk;
        %termady
        out{3} = NN1aap.*N1Npk+NNaap.*N1N1pk;
        % termNL
        out{4} = -NNpk.^3+(1-C1).*NNpk.^2+(C1).*NNpk;
        if dt_t==0 % these terms only needs to be calculated once
            % terma2_deriv
            t5 =  2*NNa.*N1Na.*cm.N1uNv+NNa.^2.*cm.lap+ 2*NNa.*NN1a.*cm.NuN1v;
        end
        % termNL_deriv
        temp = (-3*NNpk.^2+2*(1-C1).*NNpk+C1);
%         t6 = banded_dot_star(temp, cm.NuNv_flip, cm.NuNv_id);
        t6 = temp.*cm.NuNv;

        R = M_phi/tau*(out{1}-out{2}+out{3}+out{4});
        R = R*dtime-NNpk+NNp;
        dR = M_phi/tau*(t5+t6);
        dR = dR*dtime-cm.NuNv;

        %%
        % check residual and update guess
        R = R - dR*phi_initial;
        [dR, R] = kqStiffMatSetupBCID(dR,R,bcid_cp);

        dp = dR\(-R);
        phiK_cp = phiK_cp + dp;

        max_phi_R = full(max(abs(R)));
        fprintf('Phi NR Iter: %.2d -> max residual: %.2d\n',ind_check, max_phi_R);
        if (ind_check >= 100 || max(abs(R))>1e20)
            error('Phi NR method NOT converging!-Max residual: %.2d\n',max_phi_R);
        end
        ind_check = ind_check + 1;
        dt_t = dt_t+dtime;
    end

    %% Temperature (Implicit method)
    temprLHS = cm.NuNv;
    temprRHS = (cm.NuNv*tempr_cp + 3*cm.lap*tempr_cp.*dt_t + kappa*(NNpk-NNp));  %cm.NuNv\phiK
    temprRHS = temprRHS - temprLHS*tempr_initial;
    [temprLHS, temprRHS] = kqStiffMatSetupBCID(temprLHS,temprRHS,bcid_cp);

    tempr_cp_new = temprLHS\temprRHS;

    %% Theta (Implicit method)
    lap_theta = cm.lap*theta_cp;
    P2 = 10*NNpk.^3-15*NNpk.^4+6*NNpk.^5;
    P2 = P2.*s_coeff.*mag_grad_theta;

    thetaLHS = (cm.NuNv-dt_t.*M_theta.*P2.*cm.lap);
    thetaRHS = (cm.NuNv*theta_cp);
    thetaRHS = thetaRHS - thetaLHS*theta_initial;
    [thetaLHS, thetaRHS] = kqStiffMatSetupBCID(thetaLHS,thetaRHS,bcid_cp);

    theta_cp_new = thetaLHS\thetaRHS;

%     %% Tubulin concentration (Implicit method)
%     NNct = cm.NuNv*conc_t;
%     N1Nct = cm.N1uNv*conc_t;
%     NN1ct = cm.NuN1v*conc_t;
%     N2Nct = cm.N2uNv*conc_t;
%     NN2ct = cm.NuN2v*conc_t;
%     LAPct = lap*conc_t;
% 
%     sum_lap_phi = sum(LAPpk.^2);
%     NNp = cm.NuNv*phi;
%     nnpk = round(NNpk);
%     
%     term_diff = Diff.*(N1Npk.*N1Nct+NNpk.*LAPct+NN1pk.*NN1ct);
%     term_alph = alpha_t.*(N1Npk.*NNct+NNpk.*N1Nct+NN1pk.*NNct+NNpk.*NN1ct);
%     term_beta = beta_t.*NNpk.*NNct;
%     term_source = source_coeff.*LAPpk.^2./sum_lap_phi;
%     
%     conc_t_RHS = term_diff-term_alph-term_beta+term_source;
%     conc_t_RHS = (conc_t_RHS-NNct.*(NNpk-NNp)./dt_t).*dt_t./NNp+NNct;
%     conc_t_LHS = cm.NuNv;
%     bcid_t = (~nnpk);
%     [conc_t_LHS, conc_t_RHS] = StiffMatSetupBCID(conc_t_LHS, conc_t_RHS,bcid_t,zeros(lenu*lenv,1));
%     conc_t_new = conc_t_LHS\conc_t_RHS;
%     
%     %% Actin wave model
% %     actin_wave_prop

    %% iteration update
    % update variables in this iteration
    phi_cp = phiK_cp;
    theta_cp = theta_cp_new;
    tempr_cp = tempr_cp_new;
%     conct_cp = conct_cp_new;

    %% Plotting figures
    if(mod(iter,1) == 0 || iter == 1)
        phi_plot = reshape(cm.NuNv*phiK_cp,Nx,Ny);
        theta_plot = reshape(cm.NuNv*theta_cp_new,Nx,Ny);
        tempr_plot = reshape(cm.NuNv*tempr_cp_new,Nx,Ny);

        subplot(2,2,1);
        imagesc(phi_plot(2:end-1,2:end-1));
        title(sprintf('Phi plot at iter:%2d',iter));
        axis square;
        colorbar;

        subplot(2,2,2);
        displayAdaptiveGrid(ac,Coeff,Em,knotvectorU,knotvectorV,Jm,Pm,parameters,dx*nx,dy*ny);
        title(sprintf('Mesh with %d refinements',maxlev-1));
        axis square;

        subplot(2,2,3);
        imagesc(theta_plot(2:end-1,2:end-1));
        title(sprintf('Theta plot at iter:%2d',iter));
        axis square;
        colorbar;

        subplot(2,2,4);
        imagesc(tempr_plot(2:end-1,2:end-1));
        title(sprintf('T plot at iter:%2d',iter));
        axis square;
        colorbar;
% 
%         subplot(2,2,1);
%         imagesc(conct_plot);
%         title(sprintf('Tubulin plot at iter:%2d',iter));
%         axis square;
%         colorbar;

        % plot current iteration
        drawnow;

%         if(mod(iter,png_save_invl) == 0)
            try
                saveas(gcf,sprintf('./postprocessing/NeuronGrowth_%.2d.png',iter));
            catch
                fprintf('png write error skipped.\n');
            end
%         end

%         if(iter~=1 && (max(max(phi_plot(1:BC_tol,:))) > 0.5 || ...
%                 max(max(phi_plot(:,1:BC_tol))) > 0.5 || ...
%                 max(max(phi_plot(end-BC_tol:end,:))) > 0.5 || ...
%                 max(max(phi_plot(:,end-BC_tol:end))) > 0.5))
%            
%             disp('********************************************************************');
%             disp('Expanding Domain...');
% 
%             Nx = Nx+10;
%             Ny = Ny+10;
% 
%             Max_x = Max_x + 5;
%             Max_y = Max_y + 5;            
%                 
%             knotvectorU = [0,0,0,linspace(0,Nx,Nx+1),Nx,Nx,Nx].';
%             knotvectorV = [0,0,0,linspace(0,Ny,Ny+1),Ny,Ny,Ny].';
% 
%             lenu = length(knotvectorU)-2*(p-1);
%             lenv = length(knotvectorV)-2*(p-1);
% 
%             oldNuNv = cm.NuNv;
%             [cm, size_collpts] = kqCollocationDers(knotvectorU,p,knotvectorV,...
%                 q,order_deriv);
%             lap = cm.N2uNv + cm.NuN2v;
%             [lap_flip, lap_id] = extract_diags(lap);
% 
%             sz = length(lap);
% 
%             [phi,theta,theta_ori,conc_t,param,tempr,phi_initial,theta_initial,tempr_initial,bcid] ...
%                 = kqExpandDomain_Actinwave(sz,phiK,theta_new,conc_t_new,max_x,max_y,param,...
%                 tempr_new,oldNuNv,cm.NuNv);
% 
%             phiK = phi;
% 
%             toc
%             disp('********************************************************************');
%         end
    end

%     if(mod(iter,var_save_invl)==0 || iter == 0)
%         try
%             save(sprintf('./data/phi_on_cp_%2d',iter),'phi');
%             save(sprintf('./data/theta_on_cp_%2d',iter),'theta');
%             save(sprintf('./data/tempr_on_cp_%2d',iter),'tempr');
%             save(sprintf('./data/conct_on_cp_%2d',iter),'conc_t');
%         catch
%             fprintf('data write error skipped.\n');
%         end
%     end

%     phi_plot = reshape(cm.NuNv*phiK,lenu,lenv);
% 
%     if (iter < iter_stage2_begin)
%         max_x = floor(lenu/2);
%         max_y = floor(lenv/2);
% 
%     elseif (( iter>=iter_stage2_begin && iter < iter_stage3_begin) || (iter >=iter_stage45_begin) )
%             tip = sum_filter(full(phi_plot),0);
%             regionalMaxima = imregionalmax(full(tip));
%             [Max_y,Max_x] = find(regionalMaxima);
%             size_Max = length(Max_x);
%             X_dist = Max_x-lenu/2+1e-6;
%             Y_dist = Max_y-lenu/2+1e-6;
%             initial_angle = atan2(X_dist,Y_dist).';
%             
%             [theta_ori] = theta_rotate(lenu,lenv,Max_x,Max_y,initial_angle,size_Max);
%             if(mod(iter,50) == 0)
%                 subplot(2,2,3);
%                 imagesc(theta_ori);
%                 title(sprintf('Phi at iteration = %.2d',iter));
%                 axis square;
%                 colorbar;
%             end
%     elseif ( iter>=iter_stage3_begin && iter < iter_stage45_begin)
% 
%             phi_id = full(phi_plot);
%             [Nx,Ny] = size(phi_id);
%             phi_id = round(phi_id);
%             phi_sum = zeros(Nx,Ny);
% 
%             L = bwconncomp(phi_id,4);
%             S = regionprops(L,'Centroid');
%             centroids = floor(cat(1,S.Centroid));
%             
%             ID = zeros(size(phi_id));
%             dist= zeros(lenu,lenv,L.NumObjects);
% 
%             max_x = [];
%             max_y = [];
%             for k = 1:L.NumObjects
%                 ID(L.PixelIdxList{k}) = k;
%                 for i = 1:lenu
%                     for j = 1:lenv
%                         dist(i,j,k) = (ID(i,j) == k)*sqrt((i-centroids(k,1))^2+(j-centroids(k,2))^2);
%                     end
%                 end
% 
%                 dist_k = reshape(dist(:,:,k),lenu*lenv,1);
%                 [max_dist,max_index] = max(dist_k);
%                 max_x(k) = ceil(max_index/lenu);
%                 max_y(k) = rem(max_index,lenu);
% %                 if(iter == iter_stage3_begin)
% %                     x_dist = max_x(k)-centroids(k,1);
% %                     y_dist = centroids(k,2)-max_y(k);
% %                     axon_angle = atan2(x_dist,y_dist);
% %                     rotate = axon_angle;
% %                 end
% 
% %                 if ( mod(iter,rot_iter_invl) == 0 || expd_state == 1)
% %                     Rot = rand*pi/2-pi/4;
% %                     rotate_intv = Rot/rot_iter_invl;
% %                 end
% % 
% %                 theta_ori(k) = theta_rotate_guide_sector(centroids(k,1),centroids(k,2),max_x(k),max_y(k),rotate(k));
% %                 rotate(k) = rotate(k) + rotate_intv;
%             end
%             size_Max = length(max_x);
%             [theta_ori] = theta_rotate(lenu,lenv,max_x,max_y,1,size_Max);
%             
%             if(mod(iter,50) == 0)
%                 subplot(2,2,3);
%                 imagesc(theta_ori);
%                 title(sprintf('Phi at iteration = %.2d',iter));
%                 axis square;
%                 colorbar;
%             end
%     end
end

disp('All simulations complete!\n');
