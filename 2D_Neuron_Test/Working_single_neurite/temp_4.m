% IGA-collocation Implementation for 2D neuron growth
% Kuanren Qian
% 07/07/2021

%% CleanUp
close all;
clear;
clc;

%% Including Path
addpath('./IGA_collocation_algorithm');

disp('********************************************************************');
disp('2D Phase-field Neuron Growth solver using IGA-Collocation');
disp('********************************************************************');

%% Variable Initialization
% time stepping variables
dtime = 5e-3;
end_iter = 20000;

% tolerance for NR method
tol = 1e-4;

% B-spline curve order (U,V direction)
p = 3;
q = 3;
Nx = 40;
Ny = 40;

knotvectorU = [0,0,0,linspace(0,Nx,Nx+1),Nx,Nx,Nx].';
knotvectorV = [0,0,0,linspace(0,Ny,Ny+1),Ny,Ny,Ny].';
% knotvectorU = knotvectorU./max(knotvectorU);
% knotvectorV = knotvectorV./max(knotvectorV);

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
% delta = 0.2;
% epsilonb = 0.05;

% Seed size
seed_radius = 10;

% Expanding domain parameters
BC_tol = 10;
expd_coef = 1.2;

disp('Base variable - initialization done!');

%% Iterating variable initialization
% initializing phi and concentration based on neuron seed position
[phi, conc, conct] = initialize_neurite_growth(seed_radius, lenu, lenv);

% reshpae phi and concentration for calculation
phi = reshape(phi,lenu*lenv,1);
conc = reshape(conc,lenu*lenv,1);
conct = reshape(conct,lenu*lenv,1);

%% Constructing coef matrix
order_deriv = 2;    % highest order of derivatives to calculate
sprs = 1;   % sparse or not (for kqCollocationDers)
[NuNv,N1uNv,NuN1v,N1uN1v,N2uNv,NuN2v,N2uN2v,coll_p,size_collpts,Control_points,dersU] ...
    = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv,sprs);
lap = N2uNv + NuN2v;

phi = NuNv\phi;
conc = NuNv\conc;
conct = NuNv\conct;

% initializing theta and temperature
% rotate = rand*2*pi-pi;
rotate = 0;
rotate_intv = 0;
theta=rand(lenu,lenv);
for i = 1:lenu
    for j = 1:lenv
        x=i-lenu/2;
        y=j-lenv/2;
        theta(i,j) = (atan2(y,x));
        if isnan(theta(i,j))
            theta(i,j) = 1;
        end
        theta(i,j) = theta(i,j)+rotate;
        if(theta(i,j)>pi)
            theta(i,j) = theta(i,j) - 2*pi;
        elseif (theta(i,j)<-pi)
            theta(i,j) = theta(i,j) + 2*pi;
        end
        
    end
end

theta=reshape(theta,lenu*lenv,1);
theta = NuNv\theta;
tempr = zeros([lenu*lenv,1]);

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

% plotting initial phi
% set(gcf,'position',[700,100,700,900]);
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

phi = sparse(phi);
theta = sparse(theta);
tempr = sparse(tempr);
phi_ones = sparse(zeros([lenu*lenv,1]));
bcid = sparse(bcid);

dist= zeros(lenu*lenv,1);
flag_tip= zeros(lenu*lenv,1);
max_dist = 0;
max_index = 0;

% phi_actin = reshape(NuNv*phi,lenu,lenv);
% param = GetParam(phi_actin);
% actin_start = 100;

% transient iterations
for iter=1:1:end_iter
    tic
    fprintf('Progress: %.2d/%.2d\n',iter,end_iter);

    % calculating a and a*a' (aap) in the equation using theta and phi
    [a, ~, aap,pdy,pdx] = kqGetEpsilonAndAap(epsilonb,delta,phi,theta,NuNv,NuN1v,N1uNv);
    a = reshape(a,lenu*lenv,1);
    aap = reshape(aap,lenu*lenv,1);

%     teq = sparse(-1.*(NuNv*conc)+2);
    teq = 1;
    E = (alph./pix).*atan(gamma.*(teq-(NuNv*tempr)));       

    if(iter>1000)
        nnT = NuNv*theta;
        for k = 1:lenu*lenv
            if ( abs(nnT(k)) <=2.8 )
                E(k) = 0;
            end
        end
        
        subplot(2,2,3);
        imagesc(reshape(full(E),lenu,lenv)+reshape(full(NuNv*phi),lenu,lenv));
        title(sprintf('E at iteration = %.2d',iter));
        axis square;
        colorbar;
    end
    
    %% Phi (Implicit Nonlinear NR method)
    % NR method initial guess (guess current phi)
    phiK = phi;
    % residual for NR method (make sure it is larger than tolerance at the beginning)
    R = 2*tol;
    
    % magnitude of theta gradient (offset by 1e-12 to avoid division by zero)
    mag_grad_theta = sparse(sqrt((N1uNv*theta).^2+(NuN1v*theta).^2)+1e-12);
    C1 = sparse(E-0.5+6.*s_coeff.*mag_grad_theta);

    % NR method calculation
    ind_check = 0;
    NNa = NuNv*a;
    N1Na = N1uNv*a;
    NN1a = NuN1v*a;
    NNaap = NuNv*aap;
    N1Naap = N1uNv*aap;
    NN1aap = NuN1v*aap;

    dt_t = 0;
    while max(abs(R)) >= tol        
        NNpk = NuNv*phiK;
        N1Npk = N1uNv*phiK;
        NN1pk = NuN1v*phiK;   
        N1N1pk = N1uN1v*phiK;   
        LAPpk = lap*phiK;

        terma2 = 2*NNa.*N1Na.*N1Npk+NNa.^2.*LAPpk ...
            +2*NNa.*NN1a.*NN1pk;
        termadx = N1Naap.*NN1pk+NNaap.*N1N1pk;
        termady = NN1aap.*N1Npk+NNaap.*N1N1pk;     

%         termNL = -NNpk.^3+(1.5-E).*NNpk.^2+(E-0.5).*NNpk;
        termNL = -NNpk.^3+(1-C1).*NNpk.^2+(C1).*NNpk;

        terma2_deriv =  2*NNa.*N1Na.*N1uNv+NNa.^2.*lap ...
            + 2*NNa.*NN1a.*NuN1v;
        termadx_deriv = N1Naap.*NuN1v+NNaap.*N1uN1v;
        termady_deriv = NN1aap.*N1uNv+NNaap.*N1uN1v;       

%         termNL_deriv = -3*NNpk.^2+2*(1.5-E).*NNpk+(E-0.5);
%         termNL_deriv = termNL_deriv.*NuNv;
        termNL_deriv = -3*NNpk.^2+2*(1-C1).*NNpk+C1;
        termNL_deriv = termNL_deriv.*NuNv;

        R = M_phi/tau*(terma2-termadx+termady+termNL);
        R = R*dtime-NNpk+(NuNv*phi);
        dR = M_phi/tau*(terma2_deriv-termadx_deriv+termady_deriv+termNL_deriv);
        dR = dR*dtime-NuNv;

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
 
    dtime_theta = dt_t;
    dtime_conc = dt_t;
    dtime_tempr = dt_t;
    dtime_conct = dt_t;

    %% Temperature (Explicit method)
    temprLHS = NuNv;
    temprRHS = (NuNv*tempr + 3*lap*tempr.*dtime_tempr + kappa*(NuNv*phiK-NuNv*phi));
    temprRHS = temprRHS - temprLHS*tempr_initial;
    [temprLHS, temprRHS] = StiffMatSetupBCID(temprLHS, temprRHS,bcid,tempr_initial);
    tempr_new = temprLHS\temprRHS;
    
    %% Theta (Implicit method)
    lap_theta = lap*theta;
    P2 = 10*NNpk.^3-15*NNpk.^4+6*NNpk.^5;
    P2 = P2.*s_coeff.*mag_grad_theta;

    thetaLHS = (NuNv-dtime_theta.*M_theta.*P2.*lap);
    thetaRHS = (NuNv*theta);
    thetaRHS = thetaRHS - thetaLHS*theta_initial;    
    [thetaLHS, thetaRHS] = StiffMatSetupBCID(thetaLHS, thetaRHS,bcid,theta_initial);
    theta_new = thetaLHS\thetaRHS;

     %% Actin wave model
% %     if (iter >= actin_start)
%     len = length(NNpk);
%     pp = full(NNpk);
%     p = pp;
%     for i=1:len
%         if(pp(i)>0.05)
%             p(i) = 1;
%         else
%             p(i) = 0;
%         end
%     end
%     param.RegionMask = reshape(p,lenu,lenv);
%     param.phi = reshape(pp,len,1);
% 
%     % 2) set up information for ODE solver
%     T   = [0 param.dt];
%     sz  = size(param.A);
% 
%     Y_0 = [reshape(param.A,sz(1)*sz(2),1); reshape(param.H,sz(1)*sz(2),1)];
% 
%     % 3) call ODE solver
%     [~,Y] = ode45(@evolvesystem,T,Y_0,[],param);
% 
%     sz = size(param.H);
%     param.A = reshape(Y(end,1:end/2),sz);
%     param.H = reshape(Y(end,end/2+1:end),sz);
%     param.T_current = iter*param.dt;
% 
%     %% 4) update Hem
%     % 1) Change from state 0 to 1
%     H_C = max(0,(param.HTotal - sum(sum(param.H)))/param.HTotal);
%     kSA = conv2(param.H,param.G_SA,'same');
%     lmda = H_C*param.dt*kSA;
% 
%     Q = poissrnd(lmda,size(param.H));
%     len = sqrt(length(param.phi));
%     mk1 = (param.HState == 0) & Q & (param.RegionMask) & reshape(param.phi,len,len);
% 
%     % 2) Change from state 1 to 0
%     mk0 = (param.HState == 1) & (param.H < param.H_thresh_off) & (param.RegionMask);
% 
%     % 3) Update
%     % reset param values
%     param.H(mk0) = 0;
%     param.H(mk1) = param.H_thresh_on;
%     param.A(mk0) = 0;
% 
%     % reset states
%     param.HState(mk0) = 0;
%     param.HState(mk1) = 1;
% %     end

    %% iteration update
    % update variables in this iteration
    phi = phiK;
    theta = theta_new;
    tempr = tempr_new;
    
    %% Plotting figures
    if(mod(iter,50) == 0 || iter == 1)
        phi_plot = reshape(NuNv*phiK,lenu,lenv);
        subplot(2,2,1);
        imagesc(phi_plot);
        title(sprintf('Phi at iteration = %.2d',iter));
        axis square;
        colorbar;

        subplot(2,2,2);
        theta_plot = reshape(NuNv*theta_new,lenu,lenv);
        imagesc(theta_plot);
        title(sprintf('theta_plot at iteration = %.3f',rotate));
        axis square;
        colorbar;

        tempr_plot = reshape(NuNv*tempr_new,lenu,lenv);
        subplot(2,2,4);
        imagesc(tempr_plot(2:end-1,2:end-1));
        title(sprintf('Tempr at iteration = %.2d',iter));
        axis square;
        colorbar;

        % plot current iteration
        drawnow;

        if(mod(iter,100) == 0)
            saveas(gcf,sprintf('NeuronGrowth_ex1_%.2d.png',iter));
        end
        
        state = 0;
        if(iter~=1)
            if( max(max(tempr_plot(1:BC_tol,:))) > 0.4 || ...
                max(max(tempr_plot(:,1:BC_tol))) > 0.4 || ...
                max(max(tempr_plot(end-BC_tol:end,:))) > 0.4 || ...
                max(max(tempr_plot(:,end-BC_tol:end))) > 0.4)
           
                disp('********************************************************************');
                disp('Expanding Domain...');
                
                Nx = floor(Nx*expd_coef);
                Ny = floor(Ny*expd_coef);
                
                Nx = Nx+10;
                Ny = Ny+10;
                
                knotvectorU = [0,0,0,linspace(0,Nx,Nx+1),Nx,Nx,Nx].';
                knotvectorV = [0,0,0,linspace(0,Ny,Ny+1),Ny,Ny,Ny].';

                lenu = length(knotvectorU)-2*(p-1);
                lenv = length(knotvectorV)-2*(p-1);

                oldNuNv = NuNv;
                [NuNv,N1uNv,NuN1v,N1uN1v,N2uNv,NuN2v,N2uN2v,coll_p,size_collpts,Control_points] ...
                    = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv,sprs);
                lap = N2uNv + NuN2v;

                sz = length(lap);
                [phi,theta,~,tempr,phi_initial,theta_initial,~,tempr_initial,bcid] ...
                    = kqExpandDomain_4(sz,phiK,theta_new,theta_new,tempr_new,oldNuNv,NuNv,rotate);
                phiK = phi;
                
                state = 1;
                toc
                disp('********************************************************************');
            end
        end
    end
    
    rot_iter_start = 1000;
    if (iter < rot_iter_start)
        Rot = 0;
    elseif (iter>=rot_iter_start||state==1)        
        if(mod(iter,rot_iter_start)==0||iter==rot_iter_start||state==1)       
            state = 0;

            phi_full = full(reshape(NuNv*phiK,lenu,lenv));
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
            
%                     
%             if(iter == rot_iter_start)
%                 x_dist = max_x-lenu/2;
%                 y_dist = lenv/2-max_y;
%                 initial_angle = atan(x_dist/y_dist);
%             end
            
            Rot = rand*pi-pi/2;
            if abs(rotate+Rot)>=1.5
                Rot = -Rot;
            end
            
            rotate_intv = Rot/rot_iter_start;
        end
        
        theta = reshape(full(NuNv*theta),lenu,lenv);
        for i = 1:lenv
            for j = 1:lenv
                x=i-max_y-1;
                y=j-max_x;
                theta(i,j) = (atan2(y,x));
                if isnan(theta(i,j))
                    theta(i,j) = 1;
                end
                theta(i,j) = theta(i,j)+rotate;
                if(theta(i,j)>pi)
                    theta(i,j) = theta(i,j) - 2*pi;
                elseif (theta(i,j)<-pi)
                    theta(i,j) = theta(i,j) + 2*pi;
                end
            end
        end
        val = 3;
        for i = max_y-1:max_y+2
            for j = max_x-1:max_x+1
                theta(i,j) = val;
                val = -val;
            end
        end

        theta = sparse(NuNv\reshape(theta,lenu*lenv,1));
        rotate = rotate + rotate_intv;
        fprintf('Rot:%.2d, Rotate:%.4d, max_x:%.2d, max_y:%.2d\n', Rot, rotate, max_x, max_y);
    end    
    toc
end

disp('All simulations complete!\n');