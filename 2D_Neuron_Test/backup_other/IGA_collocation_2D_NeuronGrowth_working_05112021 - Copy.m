% IGA-collocation Implementation for 2D neuron growth
% Kuanren Qian
% 04/21/2021

%% CleanUp
close all;
clear;
clc;

%% Including Path
addpath('../IGA_collocation_algorithm');

disp('********************************************************************');
disp('2D Phase-field Neuron Growth solver using IGA-Collocation');
disp('********************************************************************');

%% Variable Initialization
% time stepping variables
dtime = 5e-3;
dtime_theta = 1e-3;
dtime_conc = 5e-3;
dtime_tempr = 5e-3;
dtime_conct = 5e-3;
end_iter = 20000;

multiplier = 3;
dtime = dtime*multiplier;
dtime_theta = dtime_theta*multiplier;
dtime_conc = dtime_conc*multiplier;
dtime_tempr = dtime_tempr*multiplier;
dtime_conct = dtime_conct*multiplier;

% tolerance for NR method
tol = 1e-3;

% B-spline curve order (U,V direction)
p = 3;
q = 3;
Nx = 60;
Ny = 60;
% Nx = 30;
% Ny = 30;
dx = 1;
dy = 1;
knotvectorU = [zeros([1,p]),0:Nx,ones([1,p])*Nx].';
knotvectorV = [zeros([1,p]),0:Ny,ones([1,p])*Ny].';
% setting lenu lenv this way for easier access to ghost nodes later on
lenu = length(knotvectorU)-2*(p-1);
lenv = length(knotvectorV)-2*(p-1);

% neuron growth variables
abar = 0.45;
% ab = (abar/(1+delta));
aniso = 6;
k_conc = 0.5;
kappa= 1.8;
alph = 0.9; % changing name to alph cause alpha is a function
pix=4.0*atan(1.0);
gamma = 15.0;
% gamma = 25.0;
% tau = ab;
tau = 0.3;
delta_w = 4*dx;
gamma_w = 1;
lambda = 1.0e-16;
b = 2*atanh(1-2*lambda); %88.9
W = (6*gamma_w*b)/delta_w;
M = 10;
% M_phi = (sqrt(2*W)/(6*abar))*M; %50
M_phi = 60;
M_theta = 0.5*M_phi;
D_cell = 10;
D_liquid = 1;
s_coeff = 3.e-6;

% Tubulin parameters
alpha_t = 0.001;
beta_t = 0.001;
Diff = 1e2;
source_coeff = 1.0;

% Seed size
seed_radius = 10;

disp('Base variable - initialization done!');

%% Iterating variable initialization
% initializing phi and concentration based on neuron seed position
phi = zeros([lenu,lenv]); % with ghost nodes (technically Nx by Ny, but becomes size_collpts after adding ghost points)
conc = 0.75.*ones([lenu,lenv]);
conct = zeros(lenu,lenv);
seed = (seed_radius*dx)^2;
for i=1:lenu
    for j=1:lenv
        if ((i-lenu/2)*(i-lenu/2)+(j-lenv/2)*(j-lenv/2) < seed)
            r = sqrt((i-lenu/2)*(i-lenu/2)+(j-lenv/2)*(j-lenv/2));
            phi(i,j) = 1.0;
            conc(i,j) = 0.5;
            conct(i,j) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
        end
    end
end

% reshpae phi and concentration for calculation
phi = reshape(phi,lenu*lenv,1);
conc = reshape(conc,lenu*lenv,1);
conct = reshape(conct,lenu*lenv,1);
phi_initial = reshape(phi,lenu,lenv);

%% Constructing coef matrix
order_deriv = 2;    % highest order of derivatives to calculate
sprs = 1;   % sparse or not (for kqCollocationDers)
[NuNv,N1uNv,NuN1v,N1uN1v,N2uNv,NuN2v,N2uN2v,coll_p,size_collpts,Control_points] = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv,sprs);
lap = N2uNv + NuN2v;

phi = NuNv\phi;
conc = NuNv\conc;
conct = NuNv\conct;
% initializing theta and temperature (random)
theta = NuNv\rand([lenu*lenv,1]);
tempr = zeros([lenu*lenv,1]);

% plotting initial phi
set(gcf,'position',[700,100,800,800]);
colormap parula;
% subplot(5,2,1);
% imagesc(phi_initial);
% title('Initial Phi');
% axis square;
% colorbar;
% drawnow;

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
conc = sparse(conc);
theta = sparse(theta);
tempr = sparse(tempr);
conct = sparse(conct);
phi_ones = sparse(zeros([lenu*lenv,1]));
bcid = sparse(bcid);

[phidyo,phidxo] = gradient_mat(phi_initial,Nx,Ny,dx,dy);
sq_grad_phi = sparse(reshape(phidyo.^2+phidxo.^2,lenu*lenv,1));
sum_grad_phi_ori = sum(sq_grad_phi);

delta_L = zeros(lenu*lenv,1);
term_change= ones(lenu*lenv,1);
dist= zeros(lenu*lenv,1);
flag_tip= zeros(lenu*lenv,1);
N = zeros(lenu);
% transient iterations
for iter=0:1:end_iter
    tic
    fprintf('Progress: %.2d/%.2d\n',iter,end_iter);

    % calculating a and a*a' (aap) in the equation using theta and phi
    [a, aap] = kqGetEpsilonAndAap(phi,theta,NuNv,NuN1v,N1uNv);
    a = reshape(a,lenu*lenv,1);
    aap = reshape(aap,lenu*lenv,1);

%     teq = sparse(-1.*(NuNv*conc)+1.5);
    teq = sparse(-1.*(NuNv*conc)+2);
%     if( (iter>10) )
%         
%         phi_cond = reshape(full(NuNv*phi),lenu*lenv,1);
%         for i =1:lenu*lenv
%             if(phi_cond(i)>=1e-4)
%                 x = floor(i/lenu);
%                 y = rem(i,lenu);
%                 dist(i) = sqrt((x-Nx/2)^2+(y-Ny/2)^2);
%             end
%         end
%         [maxv,maxi] = max(dist);
% 
%         conct_cond = reshape(full(NuNv*conct),lenu*lenv,1);
%         for i =1:lenu*lenv
%             x = floor(i/lenu);
%             y = rem(i,lenu);
%             if(phi_cond(i)>1e-4)
%                 if(((x-Nx/2-0.5)^2+(y-Ny/2-1.5)^2)<=(0.75*maxv)^2) % problems, need fixing
%                     flag_tip(i) = 0;
%                 else
%                     flag_tip(i) = 1;
%                 end
%             end
% 
%             %calculate delta_L
%             if( (iter>10) && (abs(phi_cond(i))>1e-4)) %&& (flag_tip(i) == 0) 
%                 r = 5e5;
%                 g = 1e-8;
% %                r = 0.1;
% %                g = -1e-8;
% 
%                 delta_L(i) = r*conct_cond(i) - g; 
%             else
%                 delta_L(i) = 1;
%             end
% 
%             if( abs(delta_L(i)) ==1 ) % just a hard-coded condition
%                 term_change(i) = 0;
%             end
%         end
%         term_change = sparse(term_change);
%         
%         E = (alph./pix).*atan(gamma.*(teq-(NuNv*tempr))).*term_change;
% 
%         % calculating teq and E term
%     %     teq = sparse(-1.*(NuNv*conc)+1.5);
%     else
%          E = (alph./pix).*atan(gamma.*(teq-(NuNv*tempr)));       
%     end
    E = (alph./pix).*atan(gamma.*(teq-(NuNv*tempr)));       

    %% Phi (Implicit Nonlinear NR method)
    % NR method initial guess (guess current phi)
    phiK = phi;
    % residual for NR method (make sure it is larger than tolerance at the beginning)
    R = 2*tol;
    
    % magnitude of theta gradient (offset by 1e-12 to avoid division by zero)
    mag_grad_theta = sparse(sqrt((N1uNv*theta).^2+(NuN1v*theta).^2)+1e-12);
    C0 = sparse(30*s_coeff*mag_grad_theta);

    % NR method calculation
    ind_check = 0;
    while max(abs(R)) >= tol        
        terma2 = 2*(NuNv*a).*(N1uNv*a).*(N1uNv*phiK)+(NuNv*a).^2.*lap*phiK ...
            +2*(NuNv*a).*(NuN1v*a).*(NuN1v*phiK);
        termadx = (N1uNv*aap).*(NuN1v*phiK)+(NuNv*aap).*(N1uN1v*phiK);
        termady = (NuN1v*aap).*(N1uNv*phiK)+(NuNv*aap).*(N1uN1v*phiK);
%         termadx = (N1uNv*a).*(NuNv*ap).*(NuN1v*phiK)+(NuNv*a).*(N1uNv*ap).*(NuN1v*phiK) ...
%             + (NuNv*a).*(NuNv*ap).*(N1uN1v*phiK);
%         termady = (NuN1v*a).*(NuNv*ap).*(N1uNv*phiK)+(NuNv*a).*(NuN1v*ap).*(N1uNv*phiK) ...
%             + (NuNv*a).*(NuNv*ap).*(N1uN1v*phiK);        
        termNL = -(NuNv*phiK).^3+(1.5-E).*(NuNv*phiK).^2+(E-0.5).*(NuNv*phiK);
        termP1theta = C0.*(NuNv*phiK).^2-2*C0.*(NuNv*phiK).^3+C0.*(NuNv*phiK).^4;

        terma2_deriv =  2*(NuNv*a).*(N1uNv*a).*N1uNv+(NuNv*a).^2.*lap ...
            + 2*(NuNv*a).*(NuN1v*a).*NuN1v;
        termadx_deriv = (N1uNv*aap).*NuN1v+(NuNv*aap).*N1uN1v;
        termady_deriv = (NuN1v*aap).*N1uNv+(NuNv*aap).*N1uN1v;
%         termadx_deriv = (N1uNv*a).*(NuNv*ap).*NuN1v+(NuNv*a).*(N1uNv*ap).*NuN1v ...
%             + (NuNv*a).*(NuNv*ap).*N1uN1v;
%         termady_deriv = (NuN1v*a).*(NuNv*ap).*N1uNv+(NuNv*a).*(NuN1v*ap).*N1uNv ...
%             + (NuNv*a).*(NuNv*ap).*N1uN1v;        
        termNL_deriv = -3*(NuNv*phiK).^2+2*(1.5-E).*(NuNv*phiK)+(E-0.5);
        termNL_deriv = termNL_deriv.*NuNv;
        termP1theta_deriv = 2*C0.*(NuNv*phiK)-6*C0.*(NuNv*phiK).^2+4*C0.*(NuNv*phiK).^3;
        termP1theta_deriv = termP1theta_deriv.*NuNv;

%         R = M_phi/tau*(terma2-termadx+termady+termNL-termP1theta);
        R = M_phi/tau*(terma2-termadx+termady+termNL);
        R = R*dtime-(NuNv*phiK)+(NuNv*phi);
        dR = M_phi/tau*(terma2_deriv-termadx_deriv+termady_deriv+termNL_deriv-termP1theta_deriv);
%         dR = M_phi/tau*(terma2_deriv-termadx_deriv+termady_deriv+termNL_deriv);
        dR = dR*dtime-NuNv;
        
        % check residual and update guess
        [dR, R] = StiffMatSetupBCID(dR, R,bcid,N);
        dp = dR\(-R);
        phiK = phiK + dp;
        
        max_phi_R = full(max(abs(R)));
%         fprintf('Phi NR Iter: %.2d -> max residual: %.2d\n',ind_check, max_phi_R);
        if (ind_check >= 100 || max(abs(R))>1e20)
            error('Phi NR method NOT converging!-Max residual: %.2d\n',max_phi_R);
        end
        ind_check = ind_check + 1;
    end
 
    %% Temperature (Implicit method)
%     temprLHS = NuNv-M_phi.*dtime_tempr.*lap;
%     temprRHS = kappa.*(NuNv*phiK-NuNv*phi)+NuNv*tempr;
%     tempr_new = temprLHS\temprRHS;
    tempr_new = NuNv\(NuNv*tempr + 5.6.*lap*tempr.*dtime_tempr + kappa*(NuNv*phiK-NuNv*phi));
%     tempr_new = NuNv\(NuNv*tempr + lap*tempr.*dtime_tempr + kappa*(NuNv*phiK-NuNv*phi));

    %% Theta (Implicit method)
    lap_theta = lap*theta;
    P2 = 10*(NuNv*phiK).^3-15*(NuNv*phiK).^4+6*(NuNv*phiK).^5;
    P2 = P2.*s_coeff.*mag_grad_theta;
%     thetaLHS = NuNv-dtime_theta.*M_theta.*P2.*lap;
%     thetaRHS = NuNv*theta;
%     [thetaLHS, thetaRHS] = StiffMatSetupBCID(thetaLHS, thetaRHS,bcid,N);
%     theta_new = thetaLHS\thetaRHS;
    theta_new = (NuNv-dtime_theta.*M_theta.*P2.*lap)\(NuNv*theta);
    
    %% Concentration (Implicit)
    term_k = (1-k_conc)./(1+(k_conc-1)*(NuNv*phiK));
    D = D_cell+(D_liquid-D_cell)*(1-(NuNv*phiK))./(1+(k_conc-1)*(NuNv*phiK));
    term_k = sparse(term_k);
    D = sparse(D);
    concL = D.*lap+D.*term_k.*(N1uNv.*(N1uNv*phiK)+NuN1v.*(NuN1v*phiK)+NuNv.*(lap*phiK));
    
%     concLHS = NuNv-concL*dtime_conc;
%     concRHS = NuNv*conc;
%     [concLHS, concRHS] = StiffMatSetupBCID(concLHS, concRHS,bcid,N);
%     conc_new = concLHS\concRHS;
    
    conc_new = (NuNv-concL*dtime_conc)\(NuNv*conc);
    
%     conc_new = NuNv\(NuNv*conc + dtime_conc.*D.*(lap*conc+term_k.*(lap*phi)));
    
    %% Tublin concentration (Implicit - incorrect implementation)
%     term_diff = (N1uNv*phiK).*N1uNv+(NuNv*phiK).*N2uNv+(NuN1v*phiK).*NuN1v+(NuNv*phiK).*NuN2v;
%     term_diff = Diff.*term_diff;
%     term_at = alpha_t.*lap;
%     term_bt = beta_t.*(NuNv*phiK).*NuNv;
%     term_e0 = source_coeff.*(sq_grad_phi./sum_grad_phi_ori);
% 
%     conctLHS = (NuNv*phiK).*NuNv-(term_diff+term_at+term_bt)*dtime_conct;
%     conctRHS = term_e0.*dtime_conct + (NuNv*phi).*(NuNv*conct);
% 
% %     conctLHS = NuNv-(term_diff-term_at-term_bt)*dtime_conct;
% %     conctRHS = term_e0.*dtime_conct + NuNv*conct;
% %     [conctLHS, conctRHS] = StiffMatSetupBCID(conctLHS, conctRHS,bcid,N);
%     conct_new = conctLHS\conctRHS;
    
    %% Plotting figures
    if(mod(iter,25) == 0 )
        phi_plot = reshape(NuNv*phiK,lenu,lenv);
        subplot(5,2,1);
        imagesc(phi_plot(2:end-1,2:end-1));
        title(sprintf('Phi at iteration = %.2d',iter));
        axis square;
        colorbar;

        subplot(5,2,2);
        E_plot = reshape(NuNv*E,lenu,lenv);
        imagesc(E_plot(2:end-1,2:end-1));
        title(sprintf('E at iteration = %.2d',iter));
        axis square;
        colorbar;

        tempr_plot = reshape(NuNv*tempr_new,lenu,lenv);
        subplot(5,2,3);
        imagesc(tempr_plot(2:end-1,2:end-1));
        title(sprintf('Tempr at iteration = %.2d',iter));
        axis square;
        colorbar;

        theta_plot = reshape(NuNv*theta_new,lenu,lenv);
        subplot(5,2,4);
        imagesc(theta_plot(2:end-1,2:end-1));
        title(sprintf('theta_plot at iteration = %.2d',iter));
        axis square;
        colorbar;

        conc_plot = reshape(NuNv*conc_new,lenu,lenv);
        subplot(5,2,5);
        imagesc(conc_plot(2:end-1,2:end-1));
        title(sprintf('Conc at iteration = %.2d',iter));
        axis square;
        colorbar;
% 
%         conct_plot = reshape(NuNv*conct_new,lenu,lenv);
%         subplot(5,2,6);
%         imagesc(conct_plot(2:end-1,2:end-1));
%         title(sprintf('conct_new at iteration = %.2d',iter));
%         axis square;
%         colorbar;

        teq_plot = reshape(teq,lenu,lenv);
        subplot(5,2,7);
        imagesc(teq_plot(2:end-1,2:end-1));
        title(sprintf('teq at %.2d',iter));
        axis square;
        colorbar;

        subplot(5,2,8);
        imagesc(reshape(NuNv*a,lenu,lenv));
        title(sprintf('a at %.2d',iter));
        axis square;
        colorbar;

        sz = sqrt(length(phi));
        fullPhi = reshape(full(NuNv*phi),sz,sz);
        [phidy,phidx]=gradient_mat_imf(fullPhi,1,1);
        pdy_plot = full(reshape(NuN1v*phiK,lenu,lenv))*20;
%         pdy_plot = reshape(phidy,lenu,lenv);
        subplot(5,2,9);
        imagesc(pdy_plot(2:end-1,2:end-1));
        title(sprintf('pdy_plot at iteration = %.2d',iter));
        axis square;
        colorbar;

        pdx_plot = full(reshape(N1uNv*phiK,lenu,lenv))*20;
        subplot(5,2,10);
        imagesc(pdx_plot(2:end-1,2:end-1));
        title(sprintf('pdx_plot at iteration = %.2d',iter));
        axis square;
        colorbar;

        t1 = NuN1v*phiK;
        t2 = N1uNv*phiK;
%         pplot = reshape(atan2(pdy_plot,pdx_plot),lenu,lenv);
        pplot = reshape(atan2(phidy,phidx),lenu,lenv);
        subplot(5,2,8);
        imagesc(pplot(2:end-1,2:end-1));
        title(sprintf('atan2_plot at iteration = %.2d',iter));
        axis square;
        colorbar;
        
        % plot current iteration
        drawnow;
    end
    
    %% iteration update
    % update variables in this iteration
    phi = phiK;
    theta = theta_new;
    conc = conc_new;
    tempr = tempr_new;
%     conct = conct_new;
% 
%     if(checkBoundaryContact(phi) == 1)
%         
%         Nx = Nx*2;
%         Ny = Ny*2;
%         knotvectorU = [zeros([1,p]),0:Nx,ones([1,p])*Nx].';
%         knotvectorV = [zeros([1,p]),0:Ny,ones([1,p])*Ny].';
%         lenu = length(knotvectorU)-2*(p-1);
%         lenv = length(knotvectorV)-2*(p-1);
% 
%         [NuNv,N1uNv,NuN1v,N1uN1v,N2uNv,NuN2v,N2uN2v,coll_p,size_collpts,Control_points] = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv,sprs);
%         lap = N2uNv + NuN2v;
%         
%         sz = length(lap);
%         [phi,theta,conc,tempr, bcid] = kqExpandDomain(sz,phi,theta,conc,tempr);
% 
%         disp('Expand');
%     end
    
    toc
end
disp('All simulations complete!\n');