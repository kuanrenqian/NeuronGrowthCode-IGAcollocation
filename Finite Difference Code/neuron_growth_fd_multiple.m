%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE-FIELD FINITE-DIFFRENCE %
%
% CODE FOR
%
% DENDRITIC SOLIDIFICATION
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%== get initial wall time:
clc
clear all
close all

time0 = clock();
format long;
%-- Simulation cell parameters:
Nx = 100;
Ny = 100;
NxNy = Nx*Ny;
dx = 1;
dy = 1;
%--- Time integration parameters:
nstep = 60000;
nprint = 100;
dtime = 1e-4;

%--- Material specific parameters:
mu = 1.0;
kappa = 1.8;
delta = 0.5;
aniso = 2.0;
abar = 0.45;
epsilonb = (abar/(1+delta));
tau = epsilonb;
M_phi = 100;
alpha = 0.9;
gamma = 10.0;
teq = 1.0;
theta0 = 0.0;
seed = 500.0;
%
pix=4.0*atan(1.0);
%--- Initialize and introduce
% initial nuclei:
[phi,tempr] = nucleus(Nx,Ny,seed);
% initializing theta and temperature
theta=rand(Nx,Ny);
theta_ori = zeros(Nx,Ny);

rot_iter_start = 1000;
rot_iter_invl = 500;

Rot = zeros(1,20);
rotate = zeros(1,20);
angle = zeros(1,20);
rotate_intv =  zeros(1,20);

rot_map = zeros(Nx,Ny);
size_Max_initial = 0;
size_Max = 0;

dist= zeros(Nx,Ny);
flag_tip= zeros(Nx,Ny);
% N = zeros(lenu*lenv,1);
max_dist = 0;
max_index = 0;

%---
%--- Evolution
%---
for istep =1:nstep
    phiold =phi;
    %---
    % calculate the laplacians
    %and epsilon:
    %---
    lap_phi = laplacian_imf(phi, dx, dy);
    %--
    lap_tempr = laplacian_imf(tempr, dx, dy);
    %--gradients of phi:
    [phidy,phidx]=gradient_mat_imf(phi,dx,dy);
    %-- calculate angle:
    atheta =atan2(phidy,phidx)+pi/2;
    %--- epsilon and its derivative:
    epsilon = epsilonb*(1.0+delta*cos(aniso*(atheta-theta)));
    epsilon_deriv = -epsilonb*aniso*delta*sin(aniso.*(atheta-theta));
    %--- first term:
    dummyx =epsilon.*epsilon_deriv.*phidx;
    [term1,~] =gradient_mat_imf(dummyx,dx,dy);
    %--- second term:
    dummyy =-epsilon.*epsilon_deriv.*phidy;
    [~,term2] =gradient_mat_imf(dummyy,dx,dy);
    %--- factor m:
    m =(alpha/pix)*atan(gamma*(teq-tempr));
    
    if(istep>=rot_iter_start)
        nnT = reshape(theta_ori,Nx*Ny,1);
        m(abs(nnT)==0) = 0;
        
        subplot(3,2,3);
        imagesc(m+phi);
        title(sprintf('E overlay with phi'));
        axis square;
        colorbar;
        
        subplot(3,2,5);
        tip = sum_filter(phi,0);
        imagesc(tip);
        title(sprintf('theta at iteration = %.2d',istep));
        axis square;
        colorbar;
    end
    
    %-- Time integration:
    phi = phi +(M_phi*dtime/tau) *(term1 +term2 + epsilon.^2 .* lap_phi + ...
        phiold.*(1.0-phiold).*(phiold - 0.5 + m));
    
    %-- evolve temperature:
    tempr = tempr + M_phi*dtime*lap_tempr + kappa*(phi-phiold);
    
    if (istep < rot_iter_start)
        max_x = floor(Nx/2);
        max_y = floor(Ny/2);
    else
        
        % neurite rotate value and get rid of center state
        if(istep>(rot_iter_start+rot_iter_invl) )
            if (mod(istep,rot_iter_invl) == 0)
                Rot = rand(1,size_Max)*2-1;
                rotate_intv = Rot/rot_iter_invl;
            end
        end
        
        % initialize max points
        if ( istep==rot_iter_start )
            tip = sum_filter(phi,0);
            tip_threshold = 1;
            while(size_Max<4)
                [Max_y,Max_x] = find(tip>tip_threshold); % arbitrary threshould
                size_Max = length(Max_x); % how many maxes
                tip_threshold = tip_threshold - 0.001;
            end
            fprintf('starting size max is : %2d', size_Max);
            X_dist = Max_x-Nx/2+1e-6;
            Y_dist = Max_y-Ny/2+1e-6;
            initial_angle = atan2(X_dist,Y_dist).';
        end
        
        % for every rot_iter_invl iterations, change max location
        %         if ( iter>rot_iter_start  && mod(iter,rot_iter_invl)==0)
        if ( istep>rot_iter_start+rot_iter_invl)
            tip = sum_filter(phi,1);
            
            regionalMaxima = imregionalmax(full(tip));
            [Max_y,Max_x] = find(regionalMaxima);
            size_Max = length(Max_x);
            X_dist = Max_x-Nx/2+1e-6;
            Y_dist = Max_y-Ny/2+1e-6;
            initial_angle = atan2(X_dist,Y_dist).';
        end
        
        [theta_ori] = theta_rotate(Nx,Ny,Max_x,Max_y,initial_angle,size_Max);
        
        subplot(3,2,6);
        imagesc(rot_map);
        title(sprintf('rot_map at iteration = %.2d',istep));
        axis square;
        colorbar;
    end
    %---- print results
    if(mod(istep,nprint) == 0 )
        fprintf('done step: %5d\n',istep);
        
        subplot(3,2,1)
        imagesc(phi)
        title("\phi")
        colorbar
        subplot(3,2,2)
        imagesc(tempr)
        title("tempr")
        colorbar
        subplot(3,2,4)
        imagesc(theta)
        title("theta")
        colorbar
        %         subplot(3,2,4)
        %         imagesc(m)
        %         title("m")
        %         colorbar
        drawnow
    end %if
    
    if(mod(istep,1000) ==0)
        str = sprintf('phi_%d.mat',istep);
        save(str, 'phi');
    end
end %istep

%--- calculate compute time:
compute_time = etime(clock(), time0);
fprintf('Compute Time: %10d\n', compute_time);