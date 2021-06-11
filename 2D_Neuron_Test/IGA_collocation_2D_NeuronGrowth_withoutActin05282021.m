% IGA-collocation Implementation for 2D neuron growth
% Kuanren Qian
% 05/27/2021

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
% dtime = 1e-3;
% dtime_theta = 1e-3;
% dtime_conc = 1e-4;
% dtime_tempr = 1e-5;
% % dtime_conct = 5e-4;
% end_iter = 200;

dtime = 5e-3;
dtime_theta = 1e-3;
dtime_conc = 5e-3;
dtime_tempr = 5e-3;
dtime_conct = 5e-3;
end_iter = 800;

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
Nx = 100;
Ny = 100;
dx = 1;
dy = 1;
knotvectorU = [zeros([1,p]),0:Nx,ones([1,p])*Nx].';
knotvectorV = [zeros([1,p]),0:Ny,ones([1,p])*Ny].';
% knotvectorU = [0,0,0,linspace(0,1,Nx),1,1,1];
% knotvectorV = [0,0,0,linspace(0,1,Ny),1,1,1];

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

% delta = 0.2;
% epsilonb = 0.005;

% delta = 0.04;
% epsilonb = 0.004;

% delta = 0.3;
% epsilonb = 0.01;
% 
% delta = 0.5;
% epsilonb = 0.008;

delta = 0.06;
epsilonb = 0.07;

% Tubulin parameters
alpha_t = 0.001;
beta_t = 0.001;
Diff = 1e2;
source_coeff = 1.0;

% Seed size
seed_radius = 15;

disp('Base variable - initialization done!');

%% Iterating variable initialization
% initializing phi and concentration based on neuron seed position
[phi, conc, conct] = initialize_neurite_growth(seed_radius, lenu, lenv, dx);

% reshpae phi and concentration for calculation
phi = reshape(phi,lenu*lenv,1);
conc = reshape(conc,lenu*lenv,1);
conct = reshape(conct,lenu*lenv,1);

%% Constructing coef matrix
order_deriv = 2;    % highest order of derivatives to calculate
sprs = 1;   % sparse or not (for kqCollocationDers)
[NuNv,N1uNv,NuN1v,N1uN1v,N2uNv,NuN2v,N2uN2v,coll_p,size_collpts,Control_points,dersU] = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv,sprs);
lap = N2uNv + NuN2v;

phi = NuNv\phi;
conc = NuNv\conc;
conct = NuNv\conct;
% initializing theta and temperature (random)
theta = NuNv\rand([lenu*lenv,1]);
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
% conct_initial  = reshape(conct,lenu*lenv,1);
theta_initial  = reshape(theta_initial,lenu*lenv,1);
tempr_initial  = reshape(tempr_initial,lenu*lenv,1);

% plotting initial phi
set(gcf,'position',[100,100,800,800]);
colormap parula;
% subplot(6,2,1);
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

phi_initial_temp = reshape(phi,lenu,lenv);
[phidyo,phidxo] = gradient_mat(phi_initial_temp,Nx,Ny,dx,dy);
sq_grad_phi = sparse(reshape(phidyo.^2+phidxo.^2,lenu*lenv,1));
sum_grad_phi_ori = sum(sq_grad_phi);

delta_L = zeros(lenu*lenv,1);
term_change= ones(lenu*lenv,1);
dist= zeros(lenu*lenv,1);
flag_tip= zeros(lenu*lenv,1);
N = zeros(lenu*lenv,1);

%% Start video writer
% writer = VideoWriter('active_wave_05272021_E_plus_A-1_0.05P'); % Define VideoWriter object
% open(writer);

% transient iterations
for iter=1:1:end_iter
    tic
    fprintf('Progress: %.2d/%.2d\n',iter,end_iter);

    % calculating a and a*a' (aap) in the equation using theta and phi
    [a, ap, aap] = kqGetEpsilonAndAap(epsilonb,delta,phi,theta,NuNv,NuN1v,N1uNv);
    a = reshape(a,lenu*lenv,1);
    aap = reshape(aap,lenu*lenv,1);

    teq = sparse(-1.*(NuNv*conc)+2);
%     teq = sparse(-1.*(NuNv*conc)+1.3);
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
    NNa = NuNv*a;
    N1Na = N1uNv*a;
    NN1a = NuN1v*a;
    
    NNap = NuNv*ap;
    N1Nap = N1uNv*ap;
    NN1ap = NuN1v*ap;
    
    NNaap = NuNv*aap;
    N1Naap = N1uNv*aap;
    NN1aap = NuN1v*aap;

    NNp = NuNv*phi;
    
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
        
%         termadx = N1Na.*NNap.*NN1pk+NNa.*N1Nap.*NN1pk ...
%             + NNa.*NNap.*N1N1pk;
%         termady = NN1a.*NNap.*N1Npk+NNa.*NN1ap.*N1Npk ...
%             + NNa.*NNap.*N1N1pk;        

        termNL = -NNpk.^3+(1.5-E).*NNpk.^2+(E-0.5).*NNpk;
        termP1theta = C0.*NNpk.^2-2*C0.*NNpk.^3+C0.*NNpk.^4;

        terma2_deriv =  2*NNa.*N1Naap.*N1uNv+NNa.^2.*lap ...
            + 2*NNa.*NN1aap.*NuN1v;
        
        termadx_deriv = N1Naap.*NuN1v+NNaap.*N1uN1v;
        termady_deriv = NN1aap.*N1uNv+NNaap.*N1uN1v;
        
%         termadx_deriv = N1Na.*NNap.*NuN1v+NNa.*N1Nap.*NuN1v ...
%             + NNa.*NNap.*N1uN1v;
%         termady_deriv = NN1a.*NNap.*N1uNv+NNa.*NN1ap.*N1uNv ...
%             + NNa.*NNap.*N1uN1v;

        termNL_deriv = -3*NNpk.^2+2*(1.5-E).*NNpk+(E-0.5);
        termNL_deriv = termNL_deriv.*NuNv;
        termP1theta_deriv = 2*C0.*NNpk-6*C0.*NNpk.^2+4*C0.*NNpk.^3;
        termP1theta_deriv = termP1theta_deriv.*NuNv;

        R = M_phi/tau*(terma2-termadx+termady+termNL-termP1theta);
        R = R*dtime-NNpk+NNp;
        dR = M_phi/tau*(terma2_deriv-termadx_deriv+termady_deriv+termNL_deriv-termP1theta_deriv);
        dR = dR*dtime-NuNv;
        
        % check residual and update guess
        R = R - dR*phi_initial;
        [dR, R] = StiffMatSetupBCID(dR, R,bcid,phi_initial);
        dp = dR\(-R);
        phiK = phiK + dp;

        max_phi_R = full(max(abs(R)));
        fprintf('Phi NR Iter: %.2d -> max residual: %.2d\n',ind_check, max_phi_R);
        if (ind_check >= 100 || max(abs(R))>1e20)
            error('Phi NR method NOT converging!-Max residual: %.2d\n',max_phi_R);
        end
        ind_check = ind_check + 1;
    end

    %% Temperature (Implicit method)
%     temprLHS = NuNv;
%     temprRHS = (NuNv*tempr + 1.3*lap*tempr.*dtime_tempr + kappa*(NuNv*phiK-NuNv*phi));
%     temprRHS = temprRHS - temprLHS*tempr_initial;
%     [temprLHS, temprRHS] = StiffMatSetupBCID(temprLHS, temprRHS,bcid,tempr_initial);
%     tempr_new = temprLHS\temprRHS;
    
    tempr_new = NuNv\(NuNv*tempr + 5.5.*lap*tempr.*dtime_tempr + kappa*(NNpk-NNp));

    %% Theta (Implicit method)
    lap_theta = lap*theta;
    P2 = 10*NNpk.^3-15*NNpk.^4+6*NNpk.^5;
    P2 = P2.*s_coeff.*mag_grad_theta;

    thetaLHS = (NuNv-dtime_theta.*M_theta.*P2.*lap);
    thetaRHS = (NuNv*theta);
    thetaRHS = thetaRHS - thetaLHS*theta_initial;
    [thetaLHS, thetaRHS] = StiffMatSetupBCID(thetaLHS, thetaRHS,bcid,theta_initial);
    theta_new = thetaLHS\thetaRHS;

%     theta_new = (NuNv-dtime_theta.*M_theta.*P2.*lap)\(NuNv*theta);

    %% Concentration (I2mplicit)
    term_k = (1-k_conc)./(1+(k_conc-1)*NNpk);
    D = D_cell+(D_liquid-D_cell)*(1-NNpk)./(1+(k_conc-1)*NNpk);
    term_k = sparse(term_k);
    D = sparse(D);
    concL = D.*lap+D.*term_k.*(N1uNv.*N1Npk+NuN1v.*NN1pk+NuNv.*LAPpk);
    
    concLHS = NuNv-concL*dtime_conc;
    concRHS = NuNv*conc;
    concRHS = concRHS - concLHS*conc_initial;
    [concLHS, concRHS] = StiffMatSetupBCID(concLHS, concRHS,bcid, conc_initial);
    conc_new = concLHS\concRHS;

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
    if(mod(iter,1)== 0 || iter == 1)
        phi_plot = reshape(NuNv*phiK,lenu,lenv);
        subplot(6,2,1);
        imagesc(phi_plot);
        title(sprintf('Phi at iteration = %.2d',iter));
        axis square;
        colorbar;

        subplot(6,2,2);
        E_plot = reshape(NuNv*E,lenu,lenv);
        imagesc(E_plot);
        title(sprintf('E at iteration = %.2d',iter));
        axis square;
        colorbar;

        tempr_plot = reshape(NuNv*tempr_new,lenu,lenv);
        subplot(6,2,3);
        imagesc(tempr_plot);
        title(sprintf('Tempr at iteration = %.2d',iter));
        axis square;
        colorbar;

        theta_plot = reshape(NuNv*theta_new,lenu,lenv);
        subplot(6,2,4);
        imagesc(theta_plot);
        title(sprintf('theta_plot at iteration = %.2d',iter));
        axis square;
        colorbar;

        conc_plot = reshape(NuNv*conc_new,lenu,lenv);
        subplot(6,2,5);
        imagesc(conc_plot);
        title(sprintf('Conc at iteration = %.2d',iter));
        axis square;
        colorbar;

        a2_plot = reshape(terma2,lenu,lenv);
        subplot(6,2,6);
        imagesc(a2_plot);
        title(sprintf('terma2 at iteration = %.2d',iter));
        axis square;
        colorbar;

        dx_plot = reshape(termadx,lenu,lenv);
        subplot(6,2,7);
        imagesc(dx_plot);
        title(sprintf('termadx at iteration = %.2d',iter));
        axis square;
        colorbar;
        
        dy_plot = reshape(termady,lenu,lenv);
        subplot(6,2,8);
        imagesc(dy_plot);
        title(sprintf('termady at iteration = %.2d',iter));
        axis square;
        colorbar;
 
        NL_plot = reshape(termNL,lenu,lenv);
        subplot(6,2,9);
        imagesc(NL_plot);
        title(sprintf('termNL at iteration = %.2d',iter));
        axis square;
        colorbar;
        
        pdx = reshape(N1uNv*phiK,lenu,lenv);
        subplot(6,2,11);
        imagesc(pdx);
        title(sprintf('pdx at iteration = %.2d',iter));
        axis square;
        colorbar;

        pdy = reshape(NuN1v*phiK,lenu,lenv);
        subplot(6,2,12);
        imagesc(pdy);
        title(sprintf('pdy at iteration = %.2d',iter));
        axis square;
        colorbar;
        
        atan2p = atan2(pdy,pdx);
        subplot(6,2,10);
        imagesc(atan2p);
        title(sprintf('atan2p at iteration = %.2d',iter));
        axis square;
        colorbar;
        
        % plot current iteration
        drawnow;
        
            
        %% Get and write the frame
%         frame = getframe(1); % Grab frame
%         writeVideo(writer,frame); % Write frame to video

    end
    
    %% iteration update
    % update variables in this iteration
    phi = phiK;
    theta = theta_new;
    conc = conc_new;
    tempr = tempr_new;
%     conct = conct_new;
    
    toc
end

% close(writer); % Close videoWriter
% fprintf('Done!\n'); % Notify user when recording is complete

disp('All simulations complete!\n');
