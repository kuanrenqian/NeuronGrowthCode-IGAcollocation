%% This script runs the iteration loop for each resolution
figure
for istep = 1:nstep
    
    %display the iteration
    str = sprintf('Iteration = %d',istep);
    disp(str);
    
    %for plotting the iteration versus N
    iteration = iteration+1;
    iter_vector = [iter_vector;iteration];
    
    if(mod(istep,250)==0)
        maxlev = parameters.maxlevel-level+1;
        
        [Dm,Pm,Em,Bvect,knotvectorU,knotvectorV,nobU,nobV,nelemU] = setBsplineGrid(maxlev,parameters,Nx,Ny,dx,dy);
        
        for multilev = 0:1:maxlev-1
            if(multilev>0)
                for j =1:bf_ct
                    bbc = bf(j,1:2);
                    bf_lev = bf(j,3);
                    EEM = Em{bf_lev,1};
                    BEM = Dm{bf_lev,1};
                    bind = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
                    supp_cells = BEM{bind,6};
                    grad = 0;
                    supp_ct = 0;
                    for i =1:size(supp_cells,2)
                        if(supp_cells(1,i)~=0)
                            supp_ct = supp_ct + 1;
                            ac_ind = EEM{supp_cells(1,i),11};
                            grad  = grad + Cell_grad(ac_ind,1);
                        end
                    end
                    
                    grad = grad/supp_ct;
                    
                    %Refinement to create next level
                    rho = parameters.rho(multilev+1);
                    if(grad>=(rho*meangrad))
                        [Dm,Em,Pm] =  Refine2D(bbc(1,1),bbc(1,2),bf_lev,Dm,Em,Pm,knotvectorU,knotvectorV,pU,pV);
                    end
                end
            end
            
            ac_ct = 0;
            bf_ct = 0;
            ac = zeros(1,2);
            bf = zeros(1,3);
            
            for lev = 1:(multilev+1)
                
                EE = Em{lev,1};
                BE = Dm{lev,1};
                sizee = size(EE,1);
                sizeb = size(BE,1);
                
                for i = 1:sizee
                    if(EE{i,4}==1)
                        ac_ct = ac_ct+1;
                        ac(ac_ct,1) = EE{i,1};
                        ac(ac_ct,2) = lev;
                        EE{i,11} = ac_ct;
                    end
                end
                
                for j = 1:sizeb,
                    if(BE{j,3}==1),
                        bf_ct = bf_ct + 1;
                        bf(bf_ct,1:2) = BE{j,1};
                        bf(bf_ct,3) = lev;
                        BE{j,10} = bf_ct;
                    end
                end
                
                Em{lev,1} = EE;
                Dm{lev,1} = BE;
            end
            
            [Jm,Coeff,Pixel] = constructAdaptiveGrid(ac,parameters,Dm,Em,X,Y,knotvectorU,knotvectorV,multilev,nobU,nobV,nelemU);
            
            Pfinal = zeros(bf_ct,2);
            for i = 1:bf_ct
                bbc = bf(i,1:2);
                bf_lev = bf(i,3);
                bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
                Pi = Pm{bf_lev,1};
                Pfinal(i,1) = Pi(bi,1);
                Pfinal(i,2) = Pi(bi,2);
            end
            
            phi_cp  = interp2(X,Y,phi,Pfinal(:,2),Pfinal(:,1));
            phi_cp(isnan(phi_cp(:))) = 0.0;
            
            cell_co = zeros(ac_ct,2);
            for i = 1:ac_ct
                cell_id = ac(i,1);
                cell_le = ac(i,2);
                EEM = Em{cell_le,1};
                cell_co(i,1) = EEM{cell_id,8};
                cell_co(i,2) = EEM{cell_id,9};
            end
            
            phi_imgg = phi.*255;
            [phidxi, phidyi] = gradient(phi_imgg);
            Idiff = sqrt(phidxi.^2+phidyi.^2);
            Cell_grad = interp2(X,Y,Idiff,cell_co(:,2),cell_co(:,1));
            meangrad = mean(Idiff(:))
        end
        
        %         if(level==1)
        %             displayAdaptiveGrid(ac,Coeff,Em,knotvectorU,knotvectorV,Jm,Pm,parameters,dx*nx,dy*ny);
        %         end
    end
    
    %Loop for calculating the gradient of phi
    px = 0;
    phi1 = zeros(Nx,Ny);
    phidx1= zeros(Nx,Ny);
    phidy1= zeros(Nx,Ny);
    
    for i = 1:Nx
        for j = 1:Ny
            px = px +1;
            ac_ind = Pixel{px,1};
            supp = Pixel{px,2};
            suppx = Pixel{px,3};
            suppy = Pixel{px,4};
            suppxx = Pixel{px,5};
            suppyy = Pixel{px,6};
            
            SB = Jm{ac_ind,1};
            ss = size(SB,1);
            f = 0;
            fx = 0;
            fy = 0;
            
            for k = 1:ss
                
                CEb = Pm{SB(k,2),1};
                BB = Dm{SB(k,2),1};
                BB_active = BB{SB(k,1),10};
                pio = CEb(SB(k,1),1);
                pj = CEb(SB(k,1),2);
                
                f = f + phi_cp(BB_active,1)*supp(k,1);
                fx = fx + phi_cp(BB_active,1)*suppx(k,1);
                fy = fy + phi_cp(BB_active,1)*suppy(k,1);
                
            end
            phidx1(j,i) = fx;
            phidy1(j,i) = fy;
        end
    end
    phidx = phidx1;
    phidy = phidy1;
    
    %store the previous iteration phi matrix
    phiold =phi;
    
    %evaluate phi*conc_tubulin
    phi_c = phi.*conc_t;
    
    %---
    % calculate the laplacians and epsilon:
    %---
    
    phi2 = reshape(phi',NxNy,1);
    
    lap_phi2 =laplacian_mat*phi2;
    
    [lap_phi]=vec2matx(lap_phi2,Nx);
    
    %--
    conc2 = reshape(conc',NxNy,1);
    
    lap_conc2 =laplacian_mat*conc2;
    
    [lap_conc]=vec2matx(lap_conc2,Nx);
    
    %--
    tempx = reshape(tempr',NxNy,1);
    
    lap_tempx =laplacian_mat*tempx;
    
    [lap_tempr]=vec2matx(lap_tempx,Nx);
    
    %--
    theta2 = reshape(theta',NxNy,1);
    
    lap_theta2 =laplacian_mat*theta2;
    
    [lap_theta]=vec2matx(lap_theta2,Nx);
    
    %--
    conct2 = reshape(conc_t',NxNy,1);
    
    lap_conct2 =laplacian_mat*conct2;
    
    [lap_conct]=vec2matx(lap_conct2,Nx);
    
    %--gradients of phi, conc, cont_tubulin and theta:
    [concdy,concdx]=gradient_mat(conc,Nx,Ny,dx,dy);
    [conctdy,conctdx]=gradient_mat(conc_t,Nx,Ny,dx,dy);
    [concvdy,concvdx]=gradient_mat(phi_c,Nx,Ny,dx,dy);
    [thetady,thetadx] = gradient_mat(theta,Nx,Ny,dx,dy);
    
    %evaluate absolute value of theta
    theta_absolute = sqrt(thetadx.^2+thetady.^2)+10-6;
    
    %evaluate d_t*phi*gradient(conc_tubulin)
    term_ct1x = phi.*conctdx;
    term_ct1y = phi.*conctdy;
    
    %calculate divergence of d_t*phi*gradient(conc_tubulin)
    [dummy,termctx] = gradient_mat(term_ct1x,Nx,Ny,dx,dy);
    [termcty,dummy] = gradient_mat(term_ct1y,Nx,Ny,dx,dy);
    
    %evaluate gradient(phi*conc_tubulin)
    phi_conc = phi.*conc_t;
    [termpcy,termpcx] = gradient_mat(phi_conc,Nx,Ny,dx,dy);
    
    %--- epsilon and its derivative:
    epsilon1 = epsilonb*(1.0+delta*cos(aniso*(pi/4)));
    epsilon_deriv1 = -epsilonb*aniso*delta*sin(aniso.*(pi/4));
    
    epsilon = zeros(Nx,Ny);
    epsilon_deriv = zeros(Nx,Ny);
    dist1 = zeros(Nx,Ny);
    radius1 = zeros(Nx,Ny);
    for i =1:Nx
        for j = 1:Ny
            
            %-- calculate angle:
            atheta =atan2(phidy(i,j),phidx(i,j));
            xtheta = 2*pi*theta(i,j);
            
            if(pi/4<=(atheta-xtheta)<=3*pi/4)
                epsilon(i,j) = epsilonb*(1.0+delta*cos(aniso*(atheta-xtheta)));
                epsilon_deriv(i,j) = -epsilonb*aniso*delta*sin(aniso.*(atheta-xtheta));
                
            elseif(-pi/4<(atheta-xtheta)<pi/4)
                epsilon(i,j) = (epsilon1/cos(pi/4)).*(cos((atheta-xtheta)));
                epsilon_deriv(i,j) = -(epsilon1/cos(pi/4)).*sin((atheta-xtheta));
                
            elseif(5*pi/4<=(atheta-xtheta)<=7*pi/4)
                epsilon(i,j) = epsilonb*(1.0+delta*cos(aniso*(atheta-xtheta)));
                epsilon_deriv(i,j) = -epsilonb*aniso*delta*sin(aniso.*(atheta-xtheta));
                
            elseif(3*pi/4<(atheta-xtheta)<5*pi/4)
                epsilon(i,j) = (epsilon1/cos(pi/4)).*(cos((atheta-xtheta)));
                epsilon_deriv(i,j) = -(epsilon1/cos(pi/4)).*sin((atheta-xtheta));
            end
        end
    end
    
    %evaluate the distance of neurites
    distance = cell(size(center,1),1);
    for kk = 1:size(center,1)
        dist1 = zeros(Nx,Ny);
        for i=1:Nx
            for j =1:Ny
                if(phi(i,j)>=0.5 && flag_center(i,j)==0.1 && kk==1)
                    dist1(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                elseif(phi(i,j)>=0.5 && flag_center(i,j)==0.2 && kk==2)
                    dist1(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                elseif(phi(i,j)>=0.5 && flag_center(i,j)==0.3 && kk==3)
                    dist1(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                elseif(phi(i,j)>=0.5 && flag_center(i,j)==0.4 && kk==4)
                    dist1(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                elseif(phi(i,j)>=0.5 && flag_center(i,j)==0.5 && kk==5)
                    dist1(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                elseif(phi(i,j)>=0.5 && flag_center(i,j)==0.6 && kk==6)
                    dist1(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                elseif(phi(i,j)>=0.5 && flag_center(i,j)==0.7 && kk==7)
                    dist1(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                elseif(phi(i,j)>=0.5 && flag_center(i,j)==0.8 && kk==8)
                    dist1(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                end
            end
        end
        distance{kk,1} = dist1;
    end
    
    %calculate longest neurite
    maxv1 = zeros(size(center,1),1);
    maxi = zeros(size(center,1),1);
    for kk = 1:size(center,1)
        dist1 = distance{kk,1};
        [maxv1(kk),maxi(kk)] = max(dist1(:));
    end
    
    %changing iteration count for selective growth
    if(iteration == iter_0)
        %flag_energy stores regions where growth is allowed
        flag_energy = zeros(Nx,Ny);
        flag_1 = zeros(Nx,Ny);
        %calculate theta where longest neurite is present
        maxv1 = zeros(size(center,1),1);
        maxi1 = zeros(size(center,1),1);
        row1 = zeros(size(center,1),1);
        col1 = zeros(size(center,1),1);
        theta_conc1 = zeros(size(center,1),1);
        
        for kk = 1:size(center,1)
            dist1 = distance{kk,1};
            [maxv1(kk),maxi1(kk)] = max(dist1(:));
            [row1(kk),col1(kk)] = ind2sub([Nx Ny],maxi1(kk));
            theta_conc1(kk,1) = atan2((row1(kk)-center(kk,1)),(col1(kk)-center(kk,2)));
            phi_before = phi;
        end
        
        
        distance2 = cell(size(center,1),1);
        for kk= 1:size(center,1)
            dist2 = zeros(Nx,Ny);
            for i =1:Nx
                for  j =1:Ny
                    %region for longest neurite
                    theta1 = atan2(i-center(kk,1),j-center(kk,2));
                    if(theta1>(theta_conc1(kk,1)-pi/8) && theta1<=(theta_conc1(kk,1)+pi/8))
                        flag_1(i,j) = 1;
                    else
                        %calculate distance for the second longest neurite
                        if(phi(i,j)>=0.5 && flag_center(i,j)==0.1 && kk==1)
                            dist2(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                        elseif(phi(i,j)>=0.5 && flag_center(i,j)==0.2 && kk==2)
                            dist2(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                        elseif(phi(i,j)>=0.5 && flag_center(i,j)==0.3 && kk==3)
                            dist2(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                        elseif(phi(i,j)>=0.5 && flag_center(i,j)==0.4 && kk==4)
                            dist2(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                        elseif(phi(i,j)>=0.5 && flag_center(i,j)==0.5 && kk==5)
                            dist2(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                        elseif(phi(i,j)>=0.5 && flag_center(i,j)==0.6 && kk==6)
                            dist2(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                        elseif(phi(i,j)>=0.5 && flag_center(i,j)==0.7 && kk==7)
                            dist2(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                        elseif(phi(i,j)>=0.5 && flag_center(i,j)==0.8 && kk==8)
                            dist2(i,j) = sqrt((i-center(kk,1))^2+(j-center(kk,2))^2);
                        end
                    end
                    distance2{kk,1} = dist2;
                end
            end
            maximum_distance = max(maxv1);
        end
        
        %second longest neurite and theta associated with it
        maxv2 = zeros(size(center,1),1);
        maxi2 = zeros(size(center,1),1);
        row2 = zeros(size(center,1),1);
        col2 = zeros(size(center,1),1);
        theta_conc2 = zeros(size(center,1),1);
        for kk = 1:size(center,1)
            dist2 = distance2{kk,1};
            [maxv2(kk),maxi2(kk)] = max(dist2(:));
            [row2(kk),col2(kk)] = ind2sub([Nx Ny],maxi2(kk));
            theta_conc2(kk,1) = atan2((row2(kk)-center(kk,1)),(col2(kk)-center(kk,2)));
            phi_before = phi;
        end
        
        %set flag_energy=1 for the longest and second longest regions
        for i=1:Nx
            for j=1:Ny
                for kk = 1:size(center,1)
                    theta1 =atan2(i-center(kk,1),j-center(kk,2));
                    if(theta1>(theta_conc1(kk,1)-pi/8) && theta1<=(theta_conc1(kk,1)+pi/8) && (sqrt((i-center(kk,1))^2+(j-center(kk,2))^2)<=3*(sqrt(seed)*dx)))
                        flag_energy(i,j) = 1;
                    end
                    
                    if(theta1>(theta_conc2(kk,1)-pi/8) && theta1<=(theta_conc2(kk,1)+pi/8)&& (sqrt((i-center(kk,1))^2+(j-center(kk,2))^2)<=3*(sqrt(seed)*dx)))
                        flag_energy(i,j) = 1;
                    end
                end
            end
        end
    end
    
    %calculate first term in phi equation:
    dummyx =epsilon.*epsilon_deriv.*phidx;
    [term1,dummy] =gradient_mat(dummyx,Nx,Ny,dx,dy);
    
    %calculate second term in phi equation:
    dummyy =-epsilon.*epsilon_deriv.*phidy;
    [dummy,term2] =gradient_mat(dummyy,Nx,Ny,dx,dy);
    
    %p(phi) calculation
    p_phi = (phi.^3).*(10-15.*phi+6.*(phi.^2));
    
    %calculate terms for theta evaluation
    dummyx = p_phi.*s_coeff.*(thetadx./theta_absolute);
    [dummy,termx] = gradient_mat(dummyx,Nx,Ny,dx,dy);
    
    dummyy = p_phi.*s_coeff.*(thetady./theta_absolute);
    [termy,dummy] = gradient_mat(dummyy,Nx,Ny,dx,dy);
    
    [termyx] = divergence(dummyx./dx,dummyy./dy);
    
    %delta_L
    delta_L = ones(Nx,Ny);
    term_change= ones(Nx,Ny);
    
    %Ematrix
    Emat = zeros(Nx,Ny);
    
    for i =1:Nx
        for  j =1:Ny
            
            %delta T
            teq = -1.*conc(i,j)+2;
            
            
            theta1 =atan2(i-center(1,1),j-center(1,2))+pi;
            if(((i-center(1,1))^2+(j-center(1,2))^2)<=(0.75*maxv1(1,1))^2)
                flag_tip = 0;
            else
                flag_tip = 1;
            end
            
            %calculate delta_L
            if(iteration>100 && iteration<iter_0 && phi(i,j)>0)
                r = 5e5;
                g = 1e-8;
                delta_L(i,j) = r*conc_t(i,j) - g;
            else
                delta_L(i,j) = 1;
            end
            
            if(delta_L(i,j)>0)
                term_change(i,j)=1;
            end
            
            %calculate E
            if(iteration>=iter_0)
                %E term in PHI evolution
                
                if(flag_energy(i,j)==0)
                    for kk =1:8
                        if(sqrt((i-center(kk,1))^2+(j-center(kk,2))^2)<=maximum_distance && kk==1 && flag_center(i,j)==0.1)
                            r = 0.1;
                            g = 0.5;
                            phi3(i,j) = 0.5;
                            delta_L(i,j) = r*conc_t(i,j) - g;
                        elseif(sqrt((i-center(kk,1))^2+(j-center(kk,2))^2)<=maximum_distance && kk==2 && flag_center(i,j)==0.2)
                            r = 0.1;
                            g = 0.5;
                            phi3(i,j) = 0.5;
                            delta_L(i,j) = r*conc_t(i,j) - g;
                        elseif(sqrt((i-center(kk,1))^2+(j-center(kk,2))^2)<=maximum_distance && kk==3 && flag_center(i,j)==0.3)
                            r = 0.1;
                            g = 0.5;
                            phi3(i,j) = 0.5;
                            delta_L(i,j) = r*conc_t(i,j) - g;
                        elseif(sqrt((i-center(kk,1))^2+(j-center(kk,2))^2)<=maximum_distance && kk==4 && flag_center(i,j)==0.4)
                            r = 0.1;
                            g = 0.5;
                            phi3(i,j) = 0.5;
                            delta_L(i,j) = r*conc_t(i,j) - g;
                        elseif(sqrt((i-center(kk,1))^2+(j-center(kk,2))^2)<=maximum_distance && kk==5 && flag_center(i,j)==0.5)
                            r = 0.1;
                            g = 0.5;
                            phi3(i,j) = 0.5;
                            delta_L(i,j) = r*conc_t(i,j) - g;
                        elseif(sqrt((i-center(kk,1))^2+(j-center(kk,2))^2)<=maximum_distance && kk==6 && flag_center(i,j)==0.6)
                            r = 0.1;
                            g = 0.5;
                            phi3(i,j) = 0.5;
                            delta_L(i,j) = r*conc_t(i,j) - g;
                        elseif(sqrt((i-center(kk,1))^2+(j-center(kk,2))^2)<=maximum_distance && kk==7 && flag_center(i,j)==0.7)
                            r = 0.1;
                            g = 0.5;
                            phi3(i,j) = 0.5;
                            delta_L(i,j) = r*conc_t(i,j) - g;
                        elseif(sqrt((i-center(kk,1))^2+(j-center(kk,2))^2)<=maximum_distance && kk==8 && flag_center(i,j)==0.8)
                            r = 0.1;
                            g = 0.5;
                            phi3(i,j) = 0.5;
                            delta_L(i,j) = r*conc_t(i,j) - g;
                        end
                    end
                elseif(flag_energy(i,j)==1 && phi(i,j)>0.5)
                    r = 5e5;
                    g = 1e-8;
                    delta_L(i,j) = r*conc_t(i,j) - g;
                end
                
                
                if(delta_L(i,j)>0)
                    term_change(i,j)=1;
                else
                    term_change(i,j)=0;
                end
            end
            
            E = (alpha/pix)*term_change(i,j)*atan(gamma*(teq-tempr(i,j)));
            Emat(i,j) = E;
            
            % calculate q(phi)
            q_phi = phi(i,j)^2*(1-phi(i,j))^2;
            
            %evolve phi:
            phi(i,j) = phi(i,j) + M_phi*(dtime/tau)*(term1(i,j) +term2(i,j) + epsilon(i,j)^2*lap_phi(i,j) + phiold(i,j)*(1.0-phiold(i,j))*(phiold(i,j) - 0.5 + E)-(30*s_coeff(i,j))*q_phi*theta_absolute(i,j));
            
            %evolve temperature:
            tempr(i,j) =tempr(i,j) + M_phi*dtime*lap_tempr(i,j) +kappa*(phi(i,j)-phiold(i,j));
            
            %evolve theta
            theta(i,j) = theta(i,j) + M_theta*(dtime_theta)*(termx(i,j)+termy(i,j));
            
            %evolve concentration
            D = D_cell+(D_liquid-D_cell)*(1-phi(i,j))/(1-phi(i,j)+k_conc*phi(i,j));
            termc = ((1-k_conc)*conc(i,j))/(1-phi(i,j)+k_conc*phi(i,j));
            
            termc_final = D*(lap_conc(i,j)+termc*lap_phi(i,j));
            conc(i,j) = conc(i,j) + (dtime_conc)*(termc_final);
            
            %evolve concentration of tubulin
            termct_diffusion = Diff.*(termctx(i,j)+termcty(i,j));
            termct_activetr = alpha_t*(termpcx(i,j)+termpcy(i,j));
            termct_degradation = beta_t*phi(i,j)*conc_t(i,j);
            termct_source = source_coeff*sq_grad_phi(i,j)/sum_grad_phi_ori;
            %termct_final = lap_conct(i,j) + lap_phi(i,j); % + alpha_t*phiold(i,j)*(conctdy(i,j)+conctdx(i,j)) + beta_t*phiold(i,j)*conc_t(i,j);
            conc_t(i,j) = conc_t(i,j)  + ((dtime_conct))*(termct_diffusion-termct_degradation+termct_source-termct_activetr);
        end
    end
    %- conc_t(i,j)*(phi(i,j)-phiold(i,j))
    phi_vec = zeros(Nx,Ny);
    phi_vec(phi>=0.5) = 1;
    
    n_vector  = [n_vector;sum(phi_vec(:))];
    
    %display images
    subplot(2,2,1)
    fig1 = imagesc(phi+phi3);
    colormap jet
    
    subplot(2,2,2)
    fig2 = imagesc(term_change);
    colormap jet
    caxis([0,1]);
    
    subplot(2,2,3)
    fig2 = imagesc(flag_energy);
    caxis([0,1]);
    colormap jet
    
    subplot(2,2,4)
    plot(iter_vector,n_vector)
    colormap jet
    drawnow
    
    %save variables
    if(mod(istep,100)==0)
        s1 = sprintf('postprocessing/phiNN4_%d.mat',iteration);
        s2 = sprintf('postprocessing/concNN4_%d.mat',iteration);
        s3 = sprintf('postprocessing/concNN4_tubulin_%d.mat',iteration);
        save(s1,'phi');
        save(s2,'conc');
        save(s3,'conc_t');
    end
    save('postprocessing/iterationNN4.mat','n_vector','iter_vector');
    
end
