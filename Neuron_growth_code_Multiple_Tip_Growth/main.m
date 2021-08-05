% IGA-collocation Implementation for 2D neuron growth
% Kuanren Qian
% 07/23/2021

%% CleanUp
close all;
clear;
%clc;

%% Including Path
addpath('./IGA_collocation_algorithm');
addpath('./Hem_algorithm');

disp('********************************************************************');
disp('2D Phase-field Neuron Growth solver using IGA-Collocation');
disp('********************************************************************');

%diary 'log_Neuron_Growth'
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
alph = 0.9;
pix=4.0*atan(1.0);
gamma = 15.0;
tau = 0.3;
M_phi = 60;
M_theta = 0.5*M_phi;
s_coeff = 0.002;

delta = 0.1;
epsilonb = 0.04;

% Seed size (radius of initial neuron seed)
seed_radius = 10;

% Expanding domain parameters
BC_tol = 10; % number of BC layers to check for expanding

rng('shuffle');

disp('Constant variable initialization done!');

%% 2D variables initialization
% Constructing coef matrix
order_deriv = 2;    % highest order of derivatives to calculate
[cm,size_collpts] = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv);
lap = cm.N2uNv + cm.NuN2v;
[lap_flip, lap_id] = extract_diags(lap);

% initializing phi and concentration based on neuron seed position
[phi, conc] = initialize_neurite_growth(seed_radius, lenu, lenv);
phi = cm.NuNv\phi;

% initializing theta and temperature
theta=cm.NuNv\reshape(rand(lenu,lenv),lenu*lenv,1);
theta_ori = zeros(lenu,lenv);
tempr = zeros(lenu*lenv,1);

% These initial variables will be used for dirichlet boundary condition
phi_initial = reshape(phi,lenu,lenv);
conc_initial = reshape(conc,lenu,lenv);
theta_initial = reshape(theta,lenu,lenv);
tempr_initial = reshape(tempr,lenu,lenv);
for i = 2:lenu-1
    for j = 2:lenv-1
        phi_initial(i,j) = 0;
        conc_initial(i,j) = 0;
        theta_initial(i,j) = 0;
        tempr_initial(i,j) = 0;
    end
end
phi_initial = reshape(phi_initial,lenu*lenv,1);
conc_initial = reshape(conc_initial,lenu*lenv,1);
theta_initial  = reshape(theta_initial,lenu*lenv,1);
tempr_initial  = reshape(tempr_initial,lenu*lenv,1);

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

% Multi-stage iteration variable setup
stage1_end = 1500;
stage2_end = 7000;
rot_iter_invl = 500;
rotate = zeros(1,20);
angle = zeros(1,20);
size_Max = 0;
Max_x = 0;
Max_y = 0;
phi_actin = reshape(cm.NuNv*phi,lenu,lenv);
param = GetParam(phi_actin,dtime);
actin_start = rot_iter_invl*2;

% save initial variable
save('./data/phi_on_cp_initial','phi');
save('./data/theta_on_cp_initial','theta');
save('./data/tempr_on_cp_initial','tempr');

% frequency to save variables/figures during simulation
var_save_invl = 1000;
png_save_invl = 100;

% converting non-sparse variable to sparse before simulation
theta = sparse(theta);
tempr = sparse(tempr);
bcid = sparse(bcid);

% Set up figure for plotting
figure(1);
% set(gcf,'position',[100,100,700,900]);
colormap parula;

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
    
    E = (alph./pix).*atan(gamma.*(1-(cm.NuNv*tempr)));
    
    if(iter>=stage1_end)
        nnT = reshape(theta_ori,lenu*lenv,1);
        E(abs(nnT)==0) = 0;
    end
    
    %% Phi (Implicit Nonlinear NR method)
    % NR method initial guess (guess current phi)
    phiK = phi;
    % residual for NR method (make sure it is larger than tolerance at the beginning)
    R = 2*tol;
    
    % magnitude of theta gradient (offset by 1e-12 to avoid division by zero)
    mag_grad_theta = sparse(sqrt((cm.N1uNv*theta).^2+(cm.NuN1v*theta).^2)+1e-12);
    C1 = sparse(E-0.5+6.*s_coeff.*mag_grad_theta);
    
    % calculate these variables in advance to save cost
    ind_check = 0;
    NNa = cm.NuNv*a;
    N1Na = cm.N1uNv*a;
    NN1a = cm.NuN1v*a;
    NNaap = cm.NuNv*aap;
    N1Naap = cm.N1uNv*aap;
    NN1aap = cm.NuN1v*aap;
    
    % out is for extracting variables from parfor
    out = cell(4,1);
    dt_t = 0;
    
    while max(abs(R)) >= tol
        NNpk = cm.NuNv*phiK;
        N1Npk = cm.N1uNv*phiK;
        NN1pk = cm.NuN1v*phiK;
        N1N1pk = cm.N1uN1v*phiK;
        LAPpk = lap*phiK;
        
        % Note: parfor seems counter-productive here as Matlab should
        % automatically use threads for most matrix and vector operations
        % like dot products
        
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
        %t6 = temp.*NuNv;
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
    
    %% Actin wave model
    %     if (iter >= actin_start)
    %         len = length(NNpk);
    %         pp = full(NNpk);
    %         for i=1:len
    %             if(pp(i)>0.05)
    %                 pp(i) = 1;
    %             else
    %                 pp(i) = 0;
    %             end
    %         end
    %         param.RegionMask = reshape(pp,lenu,lenv);
    %         param.phi = reshape(pp,len,1);
    %
    %         % 2) set up information for ODE solver
    %         T   = [0 param.dt];
    %         sz  = size(param.A);
    %
    %         Y_0 = [reshape(param.A,sz(1)*sz(2),1); reshape(param.H,sz(1)*sz(2),1)];
    %
    %         % 3) call ODE solver
    %         [~,Y] = ode45(@evolvesystem,T,Y_0,[],param);
    %
    %         sz = size(param.H);
    %         param.A = reshape(Y(end,1:end/2),sz);
    %         param.H = reshape(Y(end,end/2+1:end),sz);
    %         param.T_current = iter*param.dt;
    %
    %         %% 4) update Hem
    %         % 1) Change from state 0 to 1
    %         H_C = max(0,(param.HTotal - sum(sum(param.H)))/param.HTotal);
    %         kSA = conv2(param.H,param.G_SA,'same');
    %         lmda = H_C*param.dt*kSA;
    %
    %         Q = poissrnd(lmda,size(param.H));
    %         len = sqrt(length(param.phi));
    %         mk1 = (param.HState == 0) & Q & (param.RegionMask) & reshape(param.phi,len,len);
    %
    %         % 2) Change from state 1 to 0
    %         mk0 = (param.HState == 1) & (param.H < param.H_thresh_off) & (param.RegionMask);
    %
    %         % 3) Update
    %         % reset param values
    %         param.H(mk0) = 0;
    %         param.H(mk1) = param.H_thresh_on;
    %         param.A(mk0) = 0;
    %
    %         % reset states
    %         param.HState(mk0) = 0;
    %         param.HState(mk1) = 1;
    %     end
    
    %% iteration update
    % update variables in this iteration
    phi = phiK;
    theta = theta_new;
    tempr = tempr_new;
    
    %% Plotting figures
    if(mod(iter,50) == 0 || iter == 1)
        phi_plot = reshape(cm.NuNv*phiK,lenu,lenv);
        theta_plot = reshape(cm.NuNv*theta_new,lenu,lenv);
        tempr_plot = reshape(cm.NuNv*tempr_new,lenu,lenv);
        
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
        
        subplot(3,2,4);
        imagesc(tempr_plot(2:end-1,2:end-1));
        title(sprintf('Tempr at iteration = %.2d',iter));
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
        
        if(iter>stage1_end)
            tip = sum_filter(full(phi_plot));
            
            subplot(3,2,3);
            imagesc(reshape(E,lenu,lenv)+phi_plot);
            title(sprintf('E overlay with phi'));
            axis square;
            colorbar;
            
            subplot(3,2,5);
            imagesc(tip); hold on;
            plot(Max_x,Max_y,'o','MarkerEdgeColor','c');
            title(sprintf('theta at iteration = %.2d',iter));
            axis square;
            colorbar;
        end
        
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
            
            [phi,theta,theta_ori,param,tempr,phi_initial,theta_initial,tempr_initial,bcid] ...
                = kqExpandDomain_Actinwave(sz,phiK,theta_new,max_x,max_y,param,tempr_new,oldNuNv,cm.NuNv);
            phiK = phi;
            
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
            fprintf('Data write error skipped.\n');
        end
    end
    
    % Stage 1 (Lamellpodium initialization)
    if (iter < stage1_end)
        max_x = floor(lenu/2);
        max_y = floor(lenv/2);
        % Stage 2 (Dendrite outgrowth)
    elseif ( iter>=stage1_end && iter < stage2_end)
        tip = sum_filter(full(phi_plot));
        regionalMaxima = imregionalmax(full(tip));
        [Max_y,Max_x] = find(regionalMaxima);
        size_Max = length(Max_x);
        X_dist = Max_x-lenu/2+1e-6;
        Y_dist = Max_y-lenu/2+1e-6;
        initial_angle = atan2(X_dist,Y_dist).';
        
        theta_ori = theta_rotate(lenu,lenv,Max_x,Max_y,initial_angle,size_Max);
        
    elseif ( iter>=stage2_end)
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
        
        if(iter == stage2_end)
            x_dist = max_x-lenu/2;
            y_dist = lenv/2-max_y;
            axon_angle = atan2(x_dist,y_dist);
            rotate = axon_angle;
            
            % rotation range for axon energy sector
            axon_angle_up = axon_angle+pi/2;
            if(axon_angle_up>pi)
                axon_angle_up = axon_angle_up - 2*pi;
            elseif (axon_angle_up<-pi)
                axon_angle_up = axon_angle_up + 2*pi;
            end
            axon_angle_down = axon_angle-pi/2;
            if(axon_angle_down>pi)
                axon_angle_down = axon_angle_down - 2*pi;
            elseif (axon_angle_down<-pi)
                axon_angle_down = axon_angle_down + 2*pi;
            end
            
        end
        
        if ( mod(iter,rot_iter_invl) == 0 || expd_state == 1)
            Rot = rand*pi/2-pi/4;
            rotate_intv = Rot/rot_iter_invl;
        end
        
        theta_ori = theta_rotate(lenu,lenv,Max_x,Max_y,rotate,size_Max);
        rotate = rotate + rotate_intv;
        
        subplot(3,2,6);
        imagesc(theta_ori);
        title(sprintf('rot_map at iteration = %.2d',iter));
        axis square;
        colorbar;
    end
end

disp('All simulations complete!\n');