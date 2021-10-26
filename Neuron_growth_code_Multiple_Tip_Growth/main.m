% IGA-collocation Implementation for 2D neuron growth
% Kuanren Qian
% 08/17/2021

%% CleanUp
close all;
clear;
clc;

%% Including Path 
addpath('./IGA_collocation_algorithm');

disp('********************************************************************');
disp('2D Phase-field Neuron Growth solver using IGA-Collocation');
disp('********************************************************************');

% diary 'log_Neuron_Growth'
rngSeed = rng('shuffle');
save('./data/rngSeed','rngSeed');

%% Variable Initialization
% time stepping variables
dtime = 1e-2;
end_iter = 20000;

% tolerance for NR method
tol = 1e-4;

% B-spline curve order (U,V direction)
p = 3;
q = 3;
Nx = 60;
Ny = 60;

knotvectorU = [0,0,0,linspace(0,Nx,Nx+1),Nx,Nx,Nx].';
knotvectorV = [0,0,0,linspace(0,Ny,Ny+1),Ny,Ny,Ny].';

% setting lenu lenv this way for easier access to ghost nodes later on
lenu = length(knotvectorU)-2*(p-1);
lenv = length(knotvectorV)-2*(p-1);

% neuron growth variables
aniso = 6;
kappa= 4;
alph = 0.9; % changing name to alph cause alpha is a function
pix=4.0*atan(1.0);
alphOverPix = alph/pix;
gamma = 15.0;
tau = 0.3;
M_phi = 60;
M_theta = 0.5*M_phi;
s_coeff = 0.007;

delta = 0.1;
epsilonb = 0.04;

% Tubulin parameters
r = 1e6;
g = 0.1;
alpha_t = 0.001;
beta_t = 0.001;
Diff = 4;
source_coeff = 0.05;

% Seed size
seed_radius = 20;

% Expanding domain parameters
BC_tol = 10;
expd_coef = 1.2;

% initializing phi and concentration based on neuron seed position
[phi,conct] = initialize_neurite_growth(seed_radius, lenu, lenv);

disp('Base variable - initialization done!');

%% Constructing coef matrix
order_deriv = 2;    % highest order of derivatives to calculate
[cm,size_collpts] = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv);
lap = cm.N2uNv + cm.NuN2v;
[lap_flip, lap_id] = extract_diags(lap);

phi = cm.NuNv\phi;
conc_t = cm.NuNv\conct;

% initializing theta and temperature
theta=cm.NuNv\reshape(rand(lenu,lenv),lenu*lenv,1);
theta_ori = zeros(lenu,lenv);
tempr = zeros([lenu*lenv,1]);

conc_t_new = conct;

phi_initial = reshape(phi,lenu,lenv);
theta_initial = reshape(theta,lenu,lenv);
tempr_initial = reshape(tempr,lenu,lenv);
for i = 2:lenu-1
    for j = 2:lenv-1
        phi_initial(i,j) = 0;
        theta_initial(i,j) = 0;
        tempr_initial(i,j) = 0;
    end
end
phi_initial = reshape(phi_initial,lenu*lenv,1);
theta_initial  = reshape(theta_initial,lenu*lenv,1);
tempr_initial  = reshape(tempr_initial,lenu*lenv,1);

% plotting initial phi
set(gcf,'position',[100,100,800,400]);
colormap parula;

% ID for boundary location (suppress 4 edges)
% id = 1 means there is bc
bcid = zeros([lenu,lenv]);
for i = 1:lenu
    bcid(1,i) = 1;
    bcid(lenu,i) = 1;
    bcid(i,1) = 1;
    bcid(i,lenv) = 1;
end
bcid = reshape(bcid,lenu*lenv,1);

disp('Phi,conc,theta,tempr,bcid - initialization done!');
disp('********************************************************************');

%% Transient iteration computation
disp('Starting Neuron Growth Model transient iterations...');

theta = sparse(theta);
tempr = sparse(tempr);
bcid = sparse(bcid);

iter_stage2_begin = 500;
iter_stage3_begin = 6000;
iter_stage45_begin = 15000;
rot_iter_invl = 500;

Rot = zeros(1,20);
rotate = zeros(1,20);
rotate_intv = zeros(1,20);
size_Max = 0;

var_save_invl = 500;
png_save_invl = 100;

save('./data/phi_on_cp_initial','phi');
save('./data/theta_on_cp_initial','theta');
save('./data/tempr_on_cp_initial','tempr');

Max_x = 0;
Max_y = 0;
delta_L = 1;
term_change = 1;

tic
% transient iterations
iter = 0;
while iter <= end_iter
    iter = iter + 1;
    if(mod(iter,50) == 0)
        fprintf('Progress: %.2d/%.2d\n',iter,end_iter);
        toc
    end
    
    % calculating a and a*a' (aap) in the equation using theta and phi
    xtheta = cm.NuNv*theta;
    [a, ~, aap,~,~] = kqGetEpsilonAndAap(epsilonb,delta,phi,xtheta,cm.L_NuNv,...
        cm.U_NuNv,cm.NuN1v,cm.N1uNv);
    
    NNtempr = cm.NuNv*tempr;
    NNct = cm.NuNv*conc_t;

    if(iter<=iter_stage2_begin)
        E = alphOverPix*atan(gamma*(1-NNtempr));
    else
        nnT = reshape(theta_ori,lenu*lenv,1);
        delta_L = r*NNct - g;
        term_change = (regular_Heiviside_fun(delta_L));
        term_change(nnT==1)=1;

        E = alphOverPix*atan(gamma*term_change.*(1-NNtempr));
        E(abs(nnT)==0) = 0;
        
        if(mod(iter,50) == 0)
            subplot(2,3,5);
            phi_plot = reshape(cm.NuNv*phi,lenu,lenv);
            imagesc(reshape(E,lenu,lenv)+phi_plot);
            title(sprintf('E overlay with phi'));
            axis square;
            colorbar;
        end
    end

    %% Phi (Implicit Nonlinear NR method)
    % NR method initial guess (guess current phi)
    phiK = phi;
    % residual for NR method (make sure it is larger than tolerance at the beginning)
    R = 2*tol;
    
    % magnitude of theta gradient (offset by 1e-12 to avoid division by zero)
    mag_grad_theta = sparse(sqrt((cm.N1uNv*theta).*(cm.N1uNv*theta)+(cm.NuN1v*theta).*(cm.NuN1v*theta))+1e-12);
    C1 = sparse(E-0.5+6.*s_coeff.*mag_grad_theta);

    % NR method calculation
    ind_check = 0;
    NNa = cm.NuNv*a;
    N1Na = cm.N1uNv*a;
    NN1a = cm.NuN1v*a;
    NNaap = cm.NuNv*aap;
    N1Naap = cm.N1uNv*aap;
    NN1aap = cm.NuN1v*aap;

    out = cell(4,1);
    dt_t = 0;
    while max(abs(R)) >= tol
        NNpk = cm.NuNv*phiK;
        N1Npk = cm.N1uNv*phiK;
        NN1pk = cm.NuN1v*phiK;
        N1N1pk = cm.N1uN1v*phiK;
        LAPpk = lap*phiK;
        
        % term a2
        out{1}  = 2*NNa.*N1Na.*N1Npk+NNa.*NNa.*LAPpk ...
            +2*NNa.*NN1a.*NN1pk;
        % termadx
        out{2} = N1Naap.*NN1pk+NNaap.*N1N1pk;
        %termady
        out{3} = NN1aap.*N1Npk+NNaap.*N1N1pk;
        % termNL
        out{4} = -NNpk.*NNpk.*NNpk+(1-C1).*NNpk.*NNpk+(C1).*NNpk;
        if dt_t==0 % these terms only needs to be calculated once
            % terma2_deriv
            t5 =  banded_dot_star((2*NNa.*N1Na+N1Naap),cm.N1uNv_flip,cm.N1uNv_id) +...
                  banded_dot_star(NNa.*NNa,lap_flip,lap_id) +...
                  banded_dot_star((2*NNa.*NN1a-N1Naap),cm.NuN1v_flip,cm.NuN1v_id);
        end
        % termNL_deriv
        temp = (-3*NNpk.*NNpk+2*(1-C1).*NNpk+C1);
        t6 = banded_dot_star(temp, cm.NuNv_flip, cm.NuNv_id);
        
        R = M_phi/tau*(out{1}-out{2}+out{3}+out{4});
        R = R*dtime-NNpk+(cm.NuNv*phi);
        dR = M_phi/tau*(t5+t6);
        dR = dR*dtime-cm.NuNv;

        % check residual and update guess
        R = R - dR*phi_initial;
        [dR, R] = StiffMatSetupBCID(dR, R,bcid,phi_initial);
        dp = dR\(-R);
        phiK = phiK + dp;
        
        max_phi_R = full(max(abs(R)));
%         fprintf('Phi NR Iter: %.2d -> max residual: %.2d\n',ind_check, max_phi_R);
        if (ind_check >= 100 || max(abs(R))>1e20)
            error('Phi NR method NOT converging!-Max residual: %.2d\n',max_phi_R);
        end
        ind_check = ind_check + 1;
        dt_t = dt_t+dtime;
    end

    %% Temperature (Implicit method)
    temprLHS = cm.NuNv-3*dt_t*lap;
    temprRHS = kappa*(cm.NuNv*phiK-cm.NuNv*phi)+NNtempr;

    [temprLHS, temprRHS] = StiffMatSetupBCID(temprLHS, temprRHS,bcid,tempr_initial);
    tempr_new = temprLHS\temprRHS;

    %% Tubulin concentration (Implicit method)
    NNp = cm.NuNv*phi;
    nnpk = round(NNpk);

    if iter == 1
        initial_LAPpk = lap*phi;
        sum_lap_phi = sum(initial_LAPpk.*initial_LAPpk);
        save(sprintf('./data/initial_LAPpk_%2d',iter),'initial_LAPpk');
    end
    
    term_diff = Diff.*(banded_dot_star(N1Npk,cm.N1uNv_flip,cm.N1uNv_id) + ...
        banded_dot_star(NNpk,lap_flip,lap_id) + ...
        banded_dot_star(NN1pk,cm.NuN1v_flip,cm.NuN1v_id));
    term_alph = alpha_t.*(banded_dot_star(N1Npk,cm.NuNv_flip,cm.NuNv_id) + ...
        banded_dot_star(NNpk,cm.N1uNv_flip,cm.N1uNv_id) + ...
        banded_dot_star(NN1pk,cm.NuNv_flip,cm.NuNv_id) + ...
        banded_dot_star(NNpk,cm.NuN1v_flip,cm.NuN1v_id));
    term_beta = beta_t.*banded_dot_star(NNpk,cm.NuNv_flip,cm.NuNv_id);
    term_source = source_coeff/sum_lap_phi*initial_LAPpk.*initial_LAPpk;

    conc_t_RHS = dt_t*term_source-NNct.*(NNpk-NNp)+NNp.*NNct;
    conc_t_LHS = NNp.*cm.NuNv-dt_t.*(term_diff-term_alph-term_beta);
    bcid_t = (~nnpk);
    [conc_t_LHS, conc_t_RHS] = StiffMatSetupBCID(conc_t_LHS, conc_t_RHS,bcid_t,zeros(lenu*lenv,1));
    conc_t_new = conc_t_LHS\conc_t_RHS;
    
    %% iteration update
    % update variables in this iteration
    phi = phiK;
    tempr = tempr_new;
    conc_t = conc_t_new;

    %% Plotting figures
    if(mod(iter,50) == 0 || iter == 1)
        phi_plot = reshape(cm.NuNv*phiK,lenu,lenv);
        tempr_plot = reshape(cm.NuNv*tempr_new,lenu,lenv);
        conct_plot = reshape(cm.NuNv*conc_t_new,lenu,lenv);
        
        subplot(2,3,1);
        imagesc(phi_plot);
        title(sprintf('Phi at iteration = %.2d',iter));
        axis square;
        colorbar;

        subplot(2,3,3);
        imagesc(tempr_plot);
        title(sprintf('T at iteration = %.2d',iter));
        axis square;
        colorbar;

        subplot(2,3,4);
        imagesc(conct_plot);
        title(sprintf('Tubulin at iteration = %.2d',iter));
        axis square;
        colorbar;
        
        subplot(2,3,6);
        imagesc(reshape(initial_LAPpk,lenu,lenv));
        title(sprintf('initial_LAPpk at iteration = %.2d',iter));
        axis square;
        colorbar;

        % plot current iteration
        drawnow;

        if(mod(iter,png_save_invl) == 0)
            try
                saveas(gcf,sprintf('./data/NeuronGrowth_%.2d.png',iter));
            catch
                fprintf('png write error skipped.\n');
            end
        end
        
        if(iter~=1 && (max(max(phi_plot(1:BC_tol,:))) > 0.5 || ...
                max(max(phi_plot(:,1:BC_tol))) > 0.5 || ...
                max(max(phi_plot(end-BC_tol:end,:))) > 0.5 || ...
                max(max(phi_plot(:,end-BC_tol:end))) > 0.5))
           
            disp('********************************************************************');
            disp('Expanding Domain...');

            Nx = Nx+10;
            Ny = Ny+10;

            knotvectorU = [0,0,0,linspace(0,Nx,Nx+1),Nx,Nx,Nx].';
            knotvectorV = [0,0,0,linspace(0,Ny,Ny+1),Ny,Ny,Ny].';

            lenu = length(knotvectorU)-2*(p-1);
            lenv = length(knotvectorV)-2*(p-1);

            oldNuNv = cm.NuNv;
            [cm, size_collpts] = kqCollocationDers(knotvectorU,p,knotvectorV,...
                q,order_deriv);
            lap = cm.N2uNv + cm.NuN2v;
            [lap_flip, lap_id] = extract_diags(lap);

            sz = length(lap);

            [phi,theta,conc_t,tempr,initial_LAPpk,phi_initial,theta_initial,tempr_initial,bcid] ...
                = kqExpandDomain(sz,phiK,theta,conc_t_new,tempr_new,initial_LAPpk,oldNuNv,cm.NuNv);

            phiK = phi;

            toc
            disp('********************************************************************');
        end
    end
    
    if(mod(iter,var_save_invl)==0 || iter == 0)
        try
            save(sprintf('./data/phi_on_cp_%2d',iter),'phi');
            save(sprintf('./data/tempr_on_cp_%2d',iter),'tempr');
            save(sprintf('./data/conct_on_cp_%2d',iter),'conc_t');
        catch
            fprintf('data write error skipped.\n');
        end
    end
    
    phi_plot = reshape(cm.NuNv*phiK,lenu,lenv);

%     if (iter == 1)
%         max_x = floor(lenu/2);
%         max_y = floor(lenv/2);

    Max_x = 0;
    Max_y = 0;

    if (( iter>=iter_stage2_begin && iter < iter_stage3_begin) || (iter >=iter_stage45_begin) )
            tip = sum_filter(full(phi_plot),0);
            regionalMaxima = imregionalmax(full(tip));
            [Max_y,Max_x] = find(regionalMaxima);
            size_Max = length(Max_x);
            [theta_ori] = theta_rotate(lenu,lenv,Max_x,Max_y,size_Max);

    elseif ( iter>=iter_stage3_begin && iter < iter_stage45_begin)

            phi_id = full(phi_plot);
            [Nx,Ny] = size(phi_id);
            phi_id = round(phi_id);
            phi_sum = zeros(Nx,Ny);

            L = bwconncomp(phi_id,4);
            S = regionprops(L,'Centroid');
            centroids = floor(cat(1,S.Centroid));
            
            ID = zeros(size(phi_id));
            dist= zeros(lenu,lenv,L.NumObjects);

            max_x = [];
            max_y = [];
            for k = 1:L.NumObjects
                ID(L.PixelIdxList{k}) = k;
                for i = 1:lenu
                    for j = 1:lenv
                        dist(i,j,k) = (ID(i,j) == k)*sqrt((i-centroids(k,2))^2+(j-centroids(k,1))^2);
                    end
                end

%                 dist_k = reshape(dist(:,:,k),lenu*lenv,1);
                dist_k = dist(:,:,k);
%                 [max_dist,max_index] = max(dist_k);
                dist_k(dist_k<=0.9*max(dist_k)) = 0;
                tip = sum_filter(dist_k,0);
                regionalMaxima = imregionalmax(full(tip));
                [Max_y,Max_x] = find(regionalMaxima);

                max_x = [max_x,Max_x];
                max_y = [max_y,Max_y];

%                 for l = 1:length(max_index)
%                     max_x(end+1) = ceil(max_index(l)/lenu);
%                     max_y(end+1) = rem(max_index(l),lenu);
%                 end
%                 if(iter == iter_stage3_begin)
%                     x_dist = max_x(k)-centroids(k,1);
%                     y_dist = centroids(k,2)-max_y(k);
%                     axon_angle = atan2(x_dist,y_dist);
%                     rotate = axon_angle;
%                 end

%                 if ( mod(iter,rot_iter_invl) == 0 || expd_state == 1)
%                     Rot = rand*pi/2-pi/4;
%                     rotate_intv = Rot/rot_iter_invl;
%                 end
% 
%                 theta_ori(k) = theta_rotate_guide_sector(centroids(k,2),centroids(k,1),max_x(k),max_y(k),rotate(k));
%                 rotate(k) = rotate(k) + rotate_intv;
            end
            size_Max = length(max_x);
            [theta_ori] = theta_rotate(lenu,lenv,max_x,max_y,size_Max);
            
            if(mod(iter,50) == 0)
                subplot(2,3,3);
                imagesc(theta_ori);
                title(sprintf('Phi at iteration = %.2d',iter));
                axis square;
                colorbar;
            end
    end

end

disp('All simulations complete!\n');