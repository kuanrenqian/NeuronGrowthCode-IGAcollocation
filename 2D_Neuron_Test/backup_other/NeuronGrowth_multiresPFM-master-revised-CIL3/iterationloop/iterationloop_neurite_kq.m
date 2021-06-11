%% This script runs the iteration loop for each resolution

figure

for istep = 1:nstep
    
    tic
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
                
                for j = 1:sizeb
                    if(BE{j,3}==1)
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
            meangrad = mean(Idiff(:));
        end
        
%         if(level==1)
%             displayAdaptiveGrid(ac,Coeff,Em,knotvectorU,knotvectorV,Jm,Pm,parameters,dx*nx,dy*ny);
%         end
    end
    
%     %Loop for calculating the gradient of phi (moved to MultiResolution)
%     px = 0;
%     phi1 = zeros(Nx,Ny);
%     phidx1= zeros(Nx,Ny);
%     phidy1= zeros(Nx,Ny);

    %% Extract Dm for easy access
    % Preparing Dm variables for calculation
    % initialize Dm_extract based on Dm ({1,1} or {2,1})
    Dm1_length = length(Dm{1,1});
    Dm2_length = 0;
    Dm3_length = 0;
    if length(Dm) == 2
        Dm2_length = length(Dm{2,1});
    elseif length(Dm) == 3
        Dm2_length = length(Dm{2,1});
        Dm3_length = length(Dm{3,1});
    end
    Dm_extract = cell(Dm1_length+Dm2_length+Dm3_length,1);
    
    % asignning values from DM to Dm_extract 
    for i = 1:Dm1_length
        Dm_extract{i} = Dm{1,1}{i,10};
    end
    for i = 1:Dm2_length
        % handle empty
        if isempty(Dm{2,1}{i,10})
            Dm_extract{i+Dm1_length} = 0;
        else
            Dm_extract{i+Dm1_length} = Dm{2,1}{i,10};
        end
    end
    for i = 1:Dm3_length
        % handle empty
        if isempty(Dm{3,1}{i,10})
            Dm_extract{i+Dm1_length+Dm2_length} = 0;
        else
            Dm_extract{i+Dm1_length+Dm2_length} = Dm{3,1}{i,10};
        end
    end
    % convert to matrix for easy access
    Dm_extract = cell2mat(Dm_extract);
    
    %%  Optimization on gradient looping
%     % Pre-allocating variables (moved to MultiResolution)
%     NxNy = Nx*Ny;
%     ac_ind = zeros(NxNy,1);
%     suppx = cell(NxNy,1);
%     suppy = cell(NxNy,1);
%     Ind = cell(NxNy,1);

    % Assigning values to pre-allocated variables
    for i = 1:NxNy
        ac_ind(i,:) = Pixel{i,1};
        suppx{i} = Pixel{i,3};
        suppy{i} = Pixel{i,4};
        % Ind equation and temp variables are to calculate corresponding
        % index in Dm_extract, they are specifically related to how
        % Dm_extract variable is constructed (Complications come up when
        % refinement level increases
        temp = Jm{ac_ind(i,1),1}(:,2);
        temp = temp/2;
        temp(temp~=1)=0;
        Ind{i} = max((Jm{ac_ind(i,1),1}(:,2) - 2),0)*(Dm2_length+Dm1_length) + temp*Dm1_length + Jm{ac_ind(i,1),1}(:,1);
        
%     % Calculating phidx&phidy (moved to MultiResolution)
%     FFX = zeros(1,NxNy);
%     FFY = zeros(1,NxNy);

        BB_active_array = Dm_extract(Ind{i});
        PhiCP = phi_cp(BB_active_array);
        FFX(i) = sum(PhiCP.*suppx{i});
        FFY(i) = sum(PhiCP.*suppy{i});
    end
    phidx = reshape(FFX,Nx,Ny);
    phidy = reshape(FFY,Nx,Ny);

       %% Original code Backup
%     for i = 1:Nx
%         for j = 1:Ny
%             px = px +1;
%             ac_ind = Pixel{px,1};
%             supp = Pixel{px,2};
%             suppx = Pixel{px,3};
%             suppy = Pixel{px,4};
%             suppxx = Pixel{px,5};
%             suppyy = Pixel{px,6};
%             
%             SB = Jm{ac_ind,1};
%             ss = size(SB,1);
%             f = 0;
%             fx = 0;
%             fy = 0;
%             for k = 1:ss
%                 BB = Dm{SB(k,2),1};
%                 BB_active = BB{SB(k,1),10};
%                 %pio = CEb(SB(k,1),1);
%                 %pj = CEb(SB(k,1),2);
%                 f = f + phi_cp(BB_active,1)*supp(k,1);
%                 fx = fx + phi_cp(BB_active,1)*suppx(k,1);
%                 fy = fy + phi_cp(BB_active,1)*suppy(k,1);
%                 
%             end
%             phidx1(j,i) = fx;
%             phidy1(j,i) = fy;
%         end
%     end
%     phidx = phidx1;
%     phidy = phidy1;

    %%
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
    [~,termctx] = gradient_mat(term_ct1x,Nx,Ny,dx,dy);
    [termcty,~] = gradient_mat(term_ct1y,Nx,Ny,dx,dy);
    
    %evaluate gradient(phi*conc_tubulin)
    phi_conc = phi.*conc_t;
    [termpcy,termpcx] = gradient_mat(phi_conc,Nx,Ny,dx,dy);
    
    %--- epsilon and its derivative:
    epsilon1 = epsilonb*(1.0+delta*cos(aniso*(pi/4)));
    epsilon_deriv1 = -epsilonb*aniso*delta*sin(aniso.*(pi/4));
    
    %% Calculating Epsilon - Parallelized 
     [epsilon,epsilon_deriv,dist1] = computeEpsilon(phi,phidx,phidy,theta,epsilonb,delta,aniso,Nx,Ny,epsilon,epsilon_deriv,dist1);

    %% Epsilon code backup
%     % Loop over NxNy can be parallelized (parloop) 
%     % make it into a function (compute epsilon)
%     % output epsilon dist1
%     % %#codegen under function name
%     epsilon = zeros(Nx,Ny);
%     epsilon_deriv = zeros(Nx,Ny);
%     dist1 = zeros(Nx,Ny);
%     for i =1:Nx
%         for j = 1:Ny
%             
%             %-- calculate angle:
%             atheta =atan2(phidy(i,j),phidx(i,j));
%             xtheta = 2*pi*theta(i,j);
%             
%             if(pi/4<=(atheta-xtheta)<=3*pi/4)
%                 epsilon(i,j) = epsilonb*(1.0+delta*cos(aniso*(atheta-xtheta)));
%                 epsilon_deriv(i,j) = -epsilonb*aniso*delta*sin(aniso.*(atheta-xtheta));
%                 
%             elseif(-pi/4<(atheta-xtheta)<pi/4)
%                 epsilon(i,j) = (epsilon1/cos(pi/4)).*(cos((atheta-xtheta)));
%                 epsilon_deriv(i,j) = -(epsilon1/cos(pi/4)).*sin((atheta-xtheta));
%                 
%             elseif(5*pi/4<=(atheta-xtheta)<=7*pi/4)
%                 epsilon(i,j) = epsilonb*(1.0+delta*cos(aniso*(atheta-xtheta)));
%                 epsilon_deriv(i,j) = -epsilonb*aniso*delta*sin(aniso.*(atheta-xtheta));
%                 
%             elseif(3*pi/4<(atheta-xtheta)<5*pi/4)
%                 epsilon(i,j) = (epsilon1/cos(pi/4)).*(cos((atheta-xtheta)));
%                 epsilon_deriv(i,j) = -(epsilon1/cos(pi/4)).*sin((atheta-xtheta));
%             end
%             
%             %evaluate the distance of neurites
%             if(phi(i,j)>=0.5)
%                 dist1(i,j) = sqrt((i-Nx/2)^2+(j-Ny/2)^2);
%             end
%         end
%     end
    
    %%
    %calculate longest neurite
    [maxv1,maxi] = max(dist1(:));
    
    %calculate first term in phi equation:
    dummyx =epsilon.*epsilon_deriv.*phidx;
    [term1,~] =gradient_mat(dummyx,Nx,Ny,dx,dy);
    
    %calculate second term in phi equation:
    dummyy =-epsilon.*epsilon_deriv.*phidy;
    [~,term2] =gradient_mat(dummyy,Nx,Ny,dx,dy);
    
    %p(phi) calculation
    p_phi = (phi.^3).*(10-15.*phi+6.*(phi.^2));
    
    %calculate terms for theta evaluation
    dummyx = p_phi.*s_coeff.*(thetadx./theta_absolute);
    [~,termx] = gradient_mat(dummyx,Nx,Ny,dx,dy);
    
    dummyy = p_phi.*s_coeff.*(thetady./theta_absolute);
    [termy,~] = gradient_mat(dummyy,Nx,Ny,dx,dy);
    
    [termyx] = divergence(dummyx./dx,dummyy./dy);
    
    
    %% Main equation optimization (RHS)
%     %delta_L (moved to MultiResolution)
%     delta_L = ones(Nx,Ny);
%     term_change= ones(Nx,Ny);
%     %Could move to main
%     r = 5e5;
%     g = 1e-8;
    
    teq = -1.*conc+2;

    % is else statement needed here
    if(iteration>100)
        delta_L(phi>0) = r*conc_t(phi>0) - g;
    end
    term_change(delta_L>0) = 1;

    %E term in PHI evolution
    E = (alpha/pix)*term_change.*atan(gamma*(teq-tempr));
    
    % calculate q(phi)
    q_phi = phi.^2.*(1-phi).^2;

    %evolve phi:
    phi = phi + M_phi*(dtime/tau)*(term1 +term2 + epsilon.^2.*lap_phi + phiold.*(1.0-phiold).*(phiold - 0.5 + E)-(30*s_coeff).*q_phi.*theta_absolute);

    %evolve temperature:
    tempr =tempr + M_phi*dtime*lap_tempr +kappa*(phi-phiold);
    
    %evolve theta
    theta = theta + M_theta*(dtime_theta)*(termx+termy);

    %evolve concentration
    D = D_cell+(D_liquid-D_cell)*(1-phi)./(1-phi+k_conc*phi);
    termc = ((1-k_conc)*conc)./(1-phi+k_conc*phi);
    
    termc_final = D.*(lap_conc+termc.*lap_phi);
    conc= conc + (dtime_conc)*(termc_final);

    %evolve concentration of tubulin
    termct_diffusion = Diff.*(termctx+termcty);
    termct_activetr = alpha_t*(termpcx+termpcy);
    termct_degradation = beta_t*phi.*conc_t;
    termct_source = source_coeff*sq_grad_phi/sum_grad_phi_ori;
    conc_t = conc_t  + ((dtime_conct))*(termct_diffusion-termct_degradation+termct_source-termct_activetr);
    
    %% Main equation Backup
%     %delta_L
%     delta_L = ones(Nx,Ny);
%     term_change= ones(Nx,Ny);
%     %Ematrix
%     Emat = zeros(Nx,Ny);
%     
%     % ComputeRHS
%     for i =1:Nx
%         for  j =1:Ny
%             
%             %delta T
%             teq = -1.*conc(i,j)+2;
%             
%             theta1 =atan2(i-Nx/2,j-Ny/2)+pi;
%             if(((i-Nx/2)^2+(j-Ny/2)^2)<=(0.75*maxv1)^2)
%                 flag_tip = 0;
%             else
%                 flag_tip = 1;
%             end
%             
%             %calculate delta_L
%             if(iteration>100 && phi(i,j)>0)
%                 r = 5e5;
%                 g = 1e-8;
%                 delta_L(i,j) = r*conc_t(i,j) - g;
%             else
%                 delta_L(i,j) = 1;
%             end
%             
%             %term_change(i,j) = regular_Heiviside_fun(delta_L(i,j));
% 
%             if(delta_L(i,j)>0)
%                 term_change(i,j)=1;
%             end
%             
%             %E term in PHI evolution
%             E = (alpha/pix)*term_change(i,j)*atan(gamma*(teq-tempr(i,j)));
%             
%             % calculate q(phi)
%             q_phi = phi(i,j)^2*(1-phi(i,j))^2;
%             
%             %evolve phi:
%             phi(i,j) = phi(i,j) + M_phi*(dtime/tau)*(term1(i,j) +term2(i,j) + epsilon(i,j)^2*lap_phi(i,j) + phiold(i,j)*(1.0-phiold(i,j))*(phiold(i,j) - 0.5 + E)-(30*s_coeff(i,j))*q_phi*theta_absolute(i,j));
%             
%             %evolve temperature:
%             tempr(i,j) =tempr(i,j) + M_phi*dtime*lap_tempr(i,j) +kappa*(phi(i,j)-phiold(i,j));
%             
%             %evolve theta
%             theta(i,j) = theta(i,j) + M_theta*(dtime_theta)*(termx(i,j)+termy(i,j));
%             
%             %evolve concentration
%             D = D_cell+(D_liquid-D_cell)*(1-phi(i,j))/(1-phi(i,j)+k_conc*phi(i,j));
%             termc = ((1-k_conc)*conc(i,j))/(1-phi(i,j)+k_conc*phi(i,j));
%             
%             termc_final = D*(lap_conc(i,j)+termc*lap_phi(i,j));
%             conc(i,j) = conc(i,j) + (dtime_conc)*(termc_final);
%             
%             %evolve concentration of tubulin
%             termct_diffusion = Diff.*(termctx(i,j)+termcty(i,j));
%             termct_activetr = alpha_t*(termpcx(i,j)+termpcy(i,j));
%             termct_degradation = beta_t*phi(i,j)*conc_t(i,j);
%             termct_source = source_coeff*sq_grad_phi(i,j)/sum_grad_phi_ori;
%             %termct_final = lap_conct(i,j) + lap_phi(i,j); % + alpha_t*phiold(i,j)*(conctdy(i,j)+conctdx(i,j)) + beta_t*phiold(i,j)*conc_t(i,j);
%             conc_t(i,j) = conc_t(i,j)  + ((dtime_conct))*(termct_diffusion-termct_degradation+termct_source-termct_activetr);
%         end
%     end

    %%
    %- conc_t(i,j)*(phi(i,j)-phiold(i,j))
    phi_vec = zeros(Nx,Ny);
    phi_vec(phi>=0.5) = 1;
    
    n_vector  = [n_vector;sum(phi_vec(:))];
    
    %display images
    subplot(2,2,1)
    fig1 = imagesc(phi);
    colormap jet
    
    subplot(2,2,2)
    plot(iter_vector,n_vector)
    colormap jet
    
    subplot(2,2,3)
    fig2 = imagesc(conc);
    caxis([0,1]);
    colormap jet
    %     colorbar
    
    subplot(2,2,4)
    fig2 = imagesc(conc_t);
    %     ax = gca;
    %     ax.CLim = [0 1];
    colormap jet
    drawnow
    
    %save variables
    if(mod(istep,100)==0)
        s1 = sprintf('postprocessing/phi_nonconstraint2_%d.mat',iteration);
        s2 = sprintf('postprocessing/conc_nonconstraint2_%d.mat',iteration);
        s3 = sprintf('postprocessing/conc_tubulin_nonconstraint2_%d.mat',iteration);
        save(s1,'phi');
        save(s2,'conc');
        save(s3,'conc_t');
    end
    save('postprocessing/iteration_nonconstraint2.mat','n_vector','iter_vector');
    toc
end