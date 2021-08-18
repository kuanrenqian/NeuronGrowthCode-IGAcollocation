% IGA-collocation Implementation for 2D neuron growth
% Kuanren Qian
% 08/17/2021

%% CleanUp
close all;
clear;
clc;

%% Including Path
addpath('./IGA_collocation_algorithm');
addpath('./Hem_algorithm');

disp('********************************************************************');
disp('2D Phase-field Neuron Growth solver using IGA-Collocation');
disp('********************************************************************');

% diary 'log_Neuron_Growth'
rng('shuffle');

%% Variable Initialization
% time stepping variables
dtime = 5e-3;
end_iter = 20000;

% tolerance for NR method
tol = 1e-4;

% B-spline curve order (U,V direction)
p = 3;
q = 3;
Nx = 50;
Ny = 50;

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
gamma = 15.0;
tau = 0.3;
M_phi = 60;
M_theta = 0.5*M_phi;
s_coeff = 0.002;

delta = 0.1;
epsilonb = 0.04;

% Tubulin parameters
alpha_t = 0.001;
beta_t = 0.001;
Diff = 2;
source_coeff = 0.012;

% Seed size
seed_radius = 10;

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
set(gcf,'position',[100,-500,700,900]);
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

iter_stage2_begin = 2000;
iter_stage3_begin = 7000;
rot_iter_invl = 500;
phi_actin = reshape(cm.NuNv*phi,lenu,lenv);
param = GetParam(phi_actin,dtime);
actin_start = rot_iter_invl*2;

Rot = zeros(1,20);
rotate = zeros(1,20);
rotate_intv = zeros(1,20);
size_Max = 0;

% windowSize=5;  % Decide as per your requirements
% kernel=ones(windowSize)/windowSize^2;

var_save_invl = 1000;
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
for iter=1:1:end_iter

    if(mod(iter,50) == 0)
        fprintf('Progress: %.2d/%.2d\n',iter,end_iter);
        toc
        tic
    end
    % calculating a and a*a' (aap) in the equation using theta and phi
    xtheta = cm.NuNv*theta;
    [a, ~, aap,~,~] = kqGetEpsilonAndAap(epsilonb,delta,phi,xtheta,cm.L_NuNv,...
        cm.U_NuNv,cm.NuN1v,cm.N1uNv);
    
    if(iter<iter_stage2_begin)
        E = (alph./pix).*atan(gamma.*(1-(cm.NuNv*tempr)));
    else
        r = 5;
        g = 0.001;
        delta_L = r*reshape(conct_plot,lenu*lenv,1) - g;
        term_change = (regular_Heiviside_fun(delta_L));
        E = (alph./pix).*atan(term_change*gamma.*(1-(cm.NuNv*tempr)));

        nnT = reshape(theta_ori,lenu*lenv,1);
        E(abs(nnT)==0) = 0;
        
        subplot(3,2,3);
        phi_plot = reshape(cm.NuNv*phi,lenu,lenv);
        imagesc(reshape(E,lenu,lenv)+phi_plot);
        title(sprintf('E overlay with phi'));
        axis square;
        colorbar;
        
        subplot(3,2,5);
        tip = sum_filter(phi_plot,0);     
        imagesc(tip); hold on;
        plot(Max_x,Max_y,'o','MarkerEdgeColor','c');
        title(sprintf('theta at iteration = %.2d',iter));
        axis square;
        colorbar;
        
    end
    
    %% Phi (Implicit Nonlinear NR method)
    % NR method initial guess (guess current phi)
    phiK = phi;
    % residual for NR method (make sure it is larger than tolerance at the beginning)
    R = 2*tol;
    
    % magnitude of theta gradient (offset by 1e-12 to avoid division by zero)
    mag_grad_theta = sparse(sqrt((cm.N1uNv*theta).^2+(cm.NuN1v*theta).^2)+1e-12);
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
            t5 =  banded_dot_star((2*NNa.*N1Na+N1Naap),cm.N1uNv_flip,cm.N1uNv_id) +...
                  banded_dot_star(NNa.^2,lap_flip,lap_id) +...
                  banded_dot_star((2*NNa.*NN1a-N1Naap),cm.NuN1v_flip,cm.NuN1v_id);
        end
        % termNL_deriv
        temp = (-3*NNpk.^2+2*(1-C1).*NNpk+C1);
        t6 = banded_dot_star(temp, cm.NuNv_flip, cm.NuNv_id);
        
        R = M_phi/tau*(out{1}-out{2}+out{3}+out{4});
        R = R*dtime-NNpk+(cm.NuNv*phi);
        dR = M_phi/tau*(t5+t6);
        dR = dR*dtime-cm.NuNv;

        %%
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
    temprLHS = cm.NuNv;
    temprRHS = (cm.NuNv*tempr + 3*lap*tempr.*dt_t + kappa*(cm.NuNv*phiK-cm.NuNv*phi));
    temprRHS = temprRHS - temprLHS*tempr_initial;
    [temprLHS, temprRHS] = StiffMatSetupBCID(temprLHS, temprRHS,bcid,tempr_initial);
    tempr_new = temprLHS\temprRHS;

    %% Theta (Implicit method)
    lap_theta = lap*theta;
    P2 = 10*NNpk.^3-15*NNpk.^4+6*NNpk.^5;
    P2 = P2.*s_coeff.*mag_grad_theta;
    
    thetaLHS = (cm.NuNv-dt_t.*M_theta.*P2.*lap);
    thetaRHS = (cm.NuNv*theta);
    thetaRHS = thetaRHS - thetaLHS*theta_initial;
    [thetaLHS, thetaRHS] = StiffMatSetupBCID(thetaLHS, thetaRHS,bcid,theta_initial);
    theta_new = thetaLHS\thetaRHS;

    %% Tubulin concentration (Implicit method)
    NNct = cm.NuNv*conc_t;
    N1Nct = cm.N1uNv*conc_t;
    NN1ct = cm.NuN1v*conc_t;
    N2Nct = cm.N2uNv*conc_t;
    NN2ct = cm.NuN2v*conc_t;
    LAPct = lap*conc_t;

    sum_lap_phi = sum(LAPpk.^2);
    NNp = cm.NuNv*phi;
    nnpk = round(NNpk);
    
    term_diff = Diff.*(N1Npk.*N1Nct+NNpk.*LAPct+NN1pk.*NN1ct);
    term_alph = alpha_t.*(N1Npk.*NNct+NNpk.*N1Nct+NN1pk.*NNct+NNpk.*NN1ct);
    term_beta = beta_t.*NNpk.*NNct;
    term_source = source_coeff.*LAPpk.^2./sum_lap_phi;
    
    conc_t_RHS = term_diff-term_alph-term_beta+term_source;
    conc_t_RHS = (conc_t_RHS-NNct.*(NNpk-NNp)./dt_t).*dt_t./NNp+NNct;
    conc_t_LHS = cm.NuNv;
    bcid_t = (~nnpk);
    [conc_t_LHS, conc_t_RHS] = StiffMatSetupBCID(conc_t_LHS, conc_t_RHS,bcid_t,zeros(lenu*lenv,1));
    conc_t_new = conc_t_LHS\conc_t_RHS;
    
    %% Actin wave model
%     actin_wave_prop
    
    %% iteration update
    % update variables in this iteration
    phi = phiK;
    theta = theta_new;
    tempr = tempr_new;
    conc_t = conc_t_new;

    %% Plotting figures
    if(mod(iter,1) == 0 || iter == 1)
        phi_plot = reshape(cm.NuNv*phiK,lenu,lenv);
        theta_plot = reshape(cm.NuNv*theta_new,lenu,lenv);
        tempr_plot = reshape(cm.NuNv*tempr_new,lenu,lenv);
        conct_plot = reshape(cm.NuNv*conc_t_new,lenu,lenv);
        
        subplot(3,2,1);
        imagesc(phi_plot);
        title(sprintf('Phi at iteration = %.2d',iter));
        axis square;
        colorbar;

        subplot(3,2,2);
        imagesc(theta_plot);
        title(sprintf('theta_plot at iteration = %.2d',iter));
        axis square;
        colorbar;

%         subplot(3,2,3);
%         imagesc(tempr_plot(2:end-1,2:end-1));
%         title(sprintf('Tempr at iteration = %.2d',iter));
%         axis square;
%         colorbar;
        
        subplot(3,2,4);
        imagesc(conct_plot(2:end-1,2:end-1));
        title(sprintf('Tubulin at iteration = %.2d',iter));
        axis square;
        colorbar;
        
        %         subplot(3,2,5);
%         imagesc(param.H);
%         title(sprintf('H at iteration = %.2d',iter));
%         axis square;
%         colorbar;

%         subplot(3,2,6);
%         imagesc(param.A);
%         title(sprintf('A at iteration = %.2d',iter));
%         axis square;
%         colorbar;

        % plot current iteration
        drawnow;

        if(mod(iter,png_save_invl) == 0)
            try
                saveas(gcf,sprintf('./data/NeuronGrowth_%.2d.png',iter));
            catch
                fprintf('png write error skipped.\n');
            end
        end
        
        expd_state = 0;
        if(iter~=1 && (max(max(phi_plot(1:BC_tol,:))) > 0.5 || ...
                max(max(phi_plot(:,1:BC_tol))) > 0.5 || ...
                max(max(phi_plot(end-BC_tol:end,:))) > 0.5 || ...
                max(max(phi_plot(:,end-BC_tol:end))) > 0.5))
           
            disp('********************************************************************');
            disp('Expanding Domain...');

            Nx = Nx+10;
            Ny = Ny+10;

            Max_x = Max_x + 5;
            Max_y = Max_y + 5;            
                
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

            [phi,theta,theta_ori,conc_t,param,tempr,phi_initial,theta_initial,tempr_initial,bcid] ...
                = kqExpandDomain_Actinwave(sz,phiK,theta_new,conc_t_new,max_x,max_y,param,...
                tempr_new,oldNuNv,cm.NuNv);

            phiK = phi;

            expd_state = 1;
            toc
            disp('********************************************************************');
        end
    end
    
    if(mod(iter,var_save_invl)==0 || iter == 0)
        try
            save(sprintf('./data/phi_on_cp_%2d',iter),'phiK');
            save(sprintf('./data/theta_on_cp_%2d',iter),'theta_new');
            save(sprintf('./data/tempr_on_cp_%2d',iter),'tempr_new');
        catch
            fprintf('data write error skipped.\n');
        end
    end
    
    if (iter < iter_stage2_begin)
        max_x = floor(lenu/2);
        max_y = floor(lenv/2);

    elseif ( iter>=iter_stage2_begin && iter < iter_stage3_begin)
            tip = sum_filter(full(phi_plot),1);
            regionalMaxima = imregionalmax(full(tip));
            [Max_y,Max_x] = find(regionalMaxima);
            size_Max = length(Max_x);
            X_dist = Max_x-lenu/2+1e-6;
            Y_dist = Max_y-lenu/2+1e-6;
            initial_angle = atan2(X_dist,Y_dist).';
            
            [theta_ori] = theta_rotate(lenu,lenv,Max_x,Max_y,initial_angle,size_Max);
    elseif ( iter>=iter_stage3_begin)
                expd_state = 0;
                phi_full = full(reshape(cm.NuNv*phiK,lenu,lenv));
                dist= zeros(lenu,lenv);
                for i = 1:lenu
                    for j = 1:lenv
                        if(phi_full(i,j)>0.85)
                            dist(i,j) = sqrt((i-lenu/2)^2+(j-lenv/2)^2);
                        end
                    end
                end

                dist = reshape(dist,lenu*lenv,1);
                [max_dist,max_index] = max(dist);
                max_x = ceil(max_index/lenu);
                max_y = rem(max_index,lenu);

                if(iter == iter_stage3_begin)
                    x_dist = max_x-lenu/2;
                    y_dist = lenv/2-max_y;
                    axon_angle = atan2(x_dist,y_dist);
                    rotate = axon_angle;
                    
%                     axon_angle_up = axon_angle+pi/2;
%                     if(axon_angle_up>pi)
%                         axon_angle_up = axon_angle_up - 2*pi;
%                     elseif (axon_angle_up<-pi)
%                         axon_angle_up = axon_angle_up + 2*pi;
%                     end          
%                     axon_angle_down = axon_angle-pi/2;
%                     if(axon_angle_down>pi)
%                         axon_angle_down = axon_angle_down - 2*pi;
%                     elseif (axon_angle_down<-pi)
%                         axon_angle_down = axon_angle_down + 2*pi;
%                     end        
                    
                end
            
                if ( mod(iter,rot_iter_invl) == 0 || expd_state == 1)
                    Rot = rand*pi/2-pi/4;
%                     while ( (Rot+rotate)<axon_angle_down ||...
%                             (Rot+rotate)>axon_angle_up)
%                         Rot = rand*pi/2-pi/4;
%                     end
                    rotate_intv = Rot/rot_iter_invl;
                end

                theta_ori = theta_rotate_guide_sector(lenu,lenv,max_x,max_y,rotate);
                rotate = rotate + rotate_intv;

                subplot(3,2,6);
                imagesc(theta_ori);
                title(sprintf('rot_map at iteration = %.2d',iter));
                axis square;
                colorbar;
    end
end

disp('All simulations complete!\n');
