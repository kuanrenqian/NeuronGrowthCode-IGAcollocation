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
Nx = 200;
Ny = 200;
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
for i = 1:Nx
    for j = 1:Ny
        x=i-Nx/2;
        y=j-Ny/2;
        theta(i,j) = (atan2(y,x));
        if isnan(theta(i,j))
            theta(i,j) = 1;
        end
        theta(i,j) = theta(i,j)+0;
        if(theta(i,j)>pi)
            theta(i,j) = theta(i,j) - 2*pi;
        elseif (theta(i,j)<-pi)
            theta(i,j) = theta(i,j) + 2*pi;
        end
        
    end
end

dist= zeros(Nx,Ny);
flag_tip= zeros(Nx,Ny);
% N = zeros(lenu*lenv,1);
max_dist = 0;
max_index = 0;
rotate = 0.5;
rotate_intv = 0;
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
    if(istep>=4000)
        nnT = theta;
        %         E = E.*abs(nnT);
        for k = 1:Nx
            for l = 1:Ny
                if ( abs(nnT(k,l)) <=2.8 )
                    m(k,l) = 0;
                end
            end
        end
        %     elseif(istep>=2000)
        %         nnT = theta1;
        %         %         E = E.*abs(nnT);
        %         for k = 1:Nx
        %             for l =1:Ny
        %                 if ( abs(nnT(k,l)) <=5*pi/6)
        %                     m(k,l) = 0;
        %                 end
        %             end
        %         end
    end
    
    %-- Time integration:
    phi = phi +(M_phi*dtime/tau) *(term1 +term2 + epsilon.^2 .* lap_phi + ...
        phiold.*(1.0-phiold).*(phiold - 0.5 + m));
    
    %-- evolve temperature:
    tempr = tempr + M_phi*dtime*lap_tempr + kappa*(phi-phiold);
    
    rot_iter_start = 2000;
    if (istep < 2000)
        Rot = 0;
    elseif (istep>=2000||state==1)
        state = 0;
        if(mod(istep,rot_iter_start)==0||istep==rot_iter_start)
            Rot = rand*1-1;
            if (abs(rotate)>=1.5)
                Rot = -Rot;
            end
            rotate_intv = Rot/rot_iter_start;
            
            dist= zeros(Nx,Ny);
            for i = 1:Nx
                for j = 1:Ny
                    if(phi(i,j)>0.85)
                        dist(i,j) = sqrt((i-Nx/2-2)^2+(j-Ny/2)^2);
                    end
                end
            end
            [max_dist,max_index] = max(dist(:));
            max_x = ceil(max_index/Nx);
            max_y = rem(max_index,Nx);
        end
        for i = 1:Nx
            for j = 1:Ny
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
        rotate = rotate + rotate_intv;
        fprintf('Rot:%.2d, Rotate:%.4d, max_x:%.2d, max_y:%.2d\n', Rot, rotate, max_x, max_y);
    end
    
    if(istep==500)
        theta1 = theta;
    end
    
    %---- print results
    if(mod(istep,nprint) == 0 )
        fprintf('done step: %5d\n',istep);
        
        subplot(2,2,1)
        imagesc(phi)
        title("\phi")
        colorbar
        subplot(2,2,2)
        imagesc(tempr)
        title("tempr")
        colorbar
        subplot(2,2,3)
        imagesc(theta)
        title("theta")
        colorbar
        subplot(2,2,4)
        imagesc(m)
        title("m")
        colorbar
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