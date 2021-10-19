%% Solving laplace on THB splines
% calculating THB NuNv,N1uNv... and store in cm
cm = kqTHBders(Pixel,Jm,Pm,Dm);
[lap_flip, lap_id] = extract_diags(cm.lap);

% setting bcid (dirichlet boundary condition id)
bc_layers = 2;
[bcid,bcid_cp] = kqMakeBCID(Nx,bc_layers,X,Y,THBfinal);

% initializing phi, theta, and temperature on locally refined points
len_cp = length(THBfinal);
phi_cp = zeros(len_cp,1);
conct_cp = zeros(len_cp,1);
for i = 1:len_cp
    x = THBfinal(i,1);
    y = THBfinal(i,2);
    if (sqrt((x-0.5)^2+(y-0.5)^2) < seed_radius*dx)
        r = sqrt((x-0.5)^2+(y-0.5)^2);
        phi_cp(i) = 1;
        conct_cp(i) = (0.5+0.5*tanh(seed_radius-r)/2);
    end
end
theta_cp  = rand([len_cp,1]);
tempr_cp = zeros(len_cp,1);
theta_ori_cp = zeros(len_cp,1);

phi_initial = zeros(size(phi_cp));
theta_initial  = theta_cp;
tempr_initial  = zeros(size(tempr_cp));

disp('Phi,conc,theta,tempr,bcid - initialization done!');
disp('********************************************************************');

% variable initializations for different growth stages
kq_growth_stage_var_initialization

% just a initialization for easy later use
tempr_cp_new = tempr_cp;

update_invl = 10;
var_save_invl = 500;
png_save_invl = 100;

figure;
set(gcf,'position',[50,50,900,700]);
tic
% transient iterations
for iter=1:1:end_iter

    if(mod(iter,update_invl) ==  0 || iter == 1)
        fprintf('Progress: %.2d/%.2d\n',iter,end_iter);
        toc
    end

    % calculating a and a*a' (aap) in the equation using theta and phi
    [a_cp, ~, aap_cp,~,~] = kqGetEpsilonAndAap(epsilonb,delta,phi_cp,cm.NuNv*theta_cp,cm.NuNv,cm.NuN1v,cm.N1uNv);

    % run with simple E calculation
    tempr = cm.NuNv*tempr_cp;
    E = (alph./pix).*atan(gamma.*(1-tempr));
%     if(iter<=iter_stage2_begin)
%         E = (alph./pix).*atan(gamma.*(1-tempr));
%     else
%         conct = cm.NuNv*conct_cp;
%         r = 1e6;
%         g = 0.1;
%         delta_L = r*conct - g;
%         term_change = (regular_Heiviside_fun(delta_L));
%         term_change(theta_ori==1)=1;
%         E = (alph./pix).*atan(term_change.*gamma.*(1-tempr));
%         E(abs(theta_ori)==0) = 0;
%         
%         if(mod(iter,50) == 0)
%             subplot(2,2,2);
%             phi_plot  = griddata(THBfinal(:,1),THBfinal(:,2),phi_cp,cp_X,cp_Y);
% 
%             imagesc(E+phi_plot);
%             title(sprintf('E overlay with phi'));
%             axis square;
%             colorbar;
%         end
%     end

    %% Phi (Implicit Nonlinear NR method)
    % NR method initial guess (guess current phi)
    phiK_cp = phi_cp;
    NNp = cm.NuNv*phi_cp;

    % residual for NR method (make sure it is larger than tolerance at the beginning)
    R = 2*tol;

    % magnitude of theta gradient
    mag_grad_theta = sparse(sqrt((cm.N1uNv*theta_cp).^2+(cm.NuN1v*theta_cp).^2));
    C1 = E-0.5+6.*s_coeff.*mag_grad_theta;

    % NR method calculation
    ind_check = 0;
    NNa = cm.NuNv*a_cp;
    N1Na = cm.N1uNv*a_cp;
    NN1a = cm.NuN1v*a_cp;
    NNaap = cm.NuNv*aap_cp;
    N1Naap = cm.N1uNv*aap_cp;
    NN1aap = cm.NuN1v*aap_cp;

    out = cell(4,1);
    dt_t = 0;
    while max(abs(R)) >= tol
        NNpk = cm.NuNv*phiK_cp;
        N1Npk = cm.N1uNv*phiK_cp;
        NN1pk = cm.NuN1v*phiK_cp;
        N1N1pk = cm.N1uN1v*phiK_cp;
        LAPpk = cm.lap*phiK_cp;

        %% Non-optimized code for debugging
%         % term a2
%         out{1}  = 2*NNa.*N1Na.*N1Npk+NNa.^2.*LAPpk ...
%             +2*NNa.*NN1a.*NN1pk;
%         % termadx
%         out{2} = N1Naap.*NN1pk+NNaap.*N1N1pk;
%         %termady
%         out{3} = NN1aap.*N1Npk+NNaap.*N1N1pk;
%         % termNL
%         out{4} = -NNpk.^3+(1-C1).*NNpk.^2+(C1).*NNpk;
%         if dt_t==0 % these terms only needs to be calculated once
%             % terma2_deriv
%             t5 =  2*NNa.*N1Na.*cm.N1uNv+NNa.^2.*cm.lap+ 2*NNa.*NN1a.*cm.NuN1v;
%         end
%         % termNL_deriv
%         nl_term = (-3*NNpk.^2+2*(1-C1).*NNpk+C1);
% %         t6 = banded_dot_star(temp, cm.NuNv_flip, cm.NuNv_id);
%         t6 = nl_term.*cm.NuNv;

        %% Optimized code by Cosmin
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

        %%
        R = M_phi/tau*(out{1}-out{2}+out{3}+out{4});
        R = R*dtime-NNpk+NNp;
        dR = M_phi/tau*(t5+t6);
        dR = dR*dtime-cm.NuNv;


        % check residual and update guess
        [dR, R] = kqNGStiffMatSetupBCID(dR,R,phi_initial,bcid_cp);

        dp = dR\(-R);
        phiK_cp(bcid_cp==0) = phiK_cp(bcid_cp==0) + dp;

        max_phi_R = full(max(abs(R)));
%         fprintf('Phi NR Iter: %.2d -> max residual: %.6d\n',ind_check, max_phi_R);
        if (ind_check >= 100 || max(abs(R))>1e20)
            error('Phi NR method NOT converging!-Max residual: %.2d\n',max_phi_R);
        end
        ind_check = ind_check + 1;
        dt_t = dt_t+dtime;
    end

%     dt_tempr = dtime/5;
    dt_tempr = dtime;
    dt_theta = dtime;

    %% Temperature (Explicit method) Need to be changed to implicit
%     temprLHS = cm.NuNv;
%     temprRHS = (cm.NuNv*tempr_cp + cm.lap*tempr_cp.*dt_tempr + kappa*(NNpk-NNp));
%     [temprLHS, temprRHS] = kqNGStiffMatSetupBCID(temprLHS,temprRHS,tempr_initial,bcid_cp);
% 
%     tempr_sol = temprLHS\temprRHS;
%     tempr_cp_new(bcid_cp==0) = tempr_sol;

    %% Temperature (Implicit method)
    temprLHS = cm.NuNv-dt_tempr*cm.lap;
    temprRHS = kappa*(NNpk-NNp)+cm.NuNv*tempr_cp;
    [temprLHS, temprRHS] = kqNGStiffMatSetupBCID(temprLHS,temprRHS,tempr_initial,bcid_cp);

    tempr_sol = temprLHS\temprRHS;
    tempr_cp_new(bcid_cp==0) = tempr_sol;

    %% Theta (Implicit method)
% %     lap_theta = cm.lap*theta_cp;
%     P2 = 10*NNpk.^3-15*NNpk.^4+6*NNpk.^5;
%     P2 = P2.*s_coeff.*mag_grad_theta;
% 
%     thetaLHS = (cm.NuNv-dt_theta.*M_theta.*P2.*cm.lap);
%     thetaRHS = (cm.NuNv*theta_cp);
%     thetaRHS = thetaRHS - thetaLHS*theta_initial;
%     [thetaLHS, thetaRHS] = kqStiffMatSetupBCID(thetaLHS,thetaRHS,bcid_cp);
% 
%     theta_cp_new = thetaLHS\thetaRHS;

    % theta is de-coupled based on Ren 2018, rand theta gets the job done
    theta_cp_new = theta_cp;

    %% iteration update
    % update variables in this iteration
    phi_cp = phiK_cp;
%     theta_cp = theta_cp_new;
    tempr_cp = tempr_cp_new;
%     conct_cp = conct_cp_new;

    %% Plotting figures
    if(mod(iter,update_invl) ==  0 || iter == 1)
        % interpolating to get plottable data
        phi_plot  = griddata(THBfinal(:,1),THBfinal(:,2),cm.NuNv*phi_cp,X,Y);
        theta_plot  = griddata(THBfinal(:,1),THBfinal(:,2),cm.NuNv*theta_cp_new,X,Y);
        tempr_plot  = griddata(THBfinal(:,1),THBfinal(:,2),cm.NuNv*tempr_cp_new,X,Y);
        E_plot  = griddata(THBfinal(:,1),THBfinal(:,2),E,X,Y);
        phi_diff  = griddata(THBfinal(:,1),THBfinal(:,2),(NNpk-NNp),X,Y);

        subplot(2,3,1);
        imagesc(phi_plot);
        title(sprintf('Phi plot at iter:%2d',iter));
        axis square;
        colorbar;

        subplot(2,3,2);
        imagesc(phi_plot);
        title(sprintf('Phi plot at iter:%2d',iter));
        axis square;
        colorbar;
        hold on;
        displayAdaptiveGrid(ac,Coeff,Em,knotvectorU,knotvectorV,Jm,Pm,parameters,Nx,0.5);
        title(sprintf('Mesh with %d refinements',maxlev-1));
        axis square;

        subplot(2,3,3);
        imagesc(theta_plot);
        title(sprintf('Theta plot at iter:%2d',iter));
        axis square;
        colorbar;

        subplot(2,3,4);
        imagesc(tempr_plot);
        title(sprintf('T plot at iter:%2d',iter));
        axis square;
        colorbar;
        
        subplot(2,3,5);
        imagesc(E_plot);
        title(sprintf('E at iter:%2d',iter));
        axis square;
        colorbar;

        subplot(2,3,6);
        imagesc(phi_diff);
        title('Phi diff');
        axis square;
        colorbar;

        % plot current iteration
        drawnow;

        % save plot as png
        try
            saveas(gcf,sprintf('./postprocessing/NeuronGrowth_%.2d.png',iter));
        catch
            fprintf('png write error skipped.\n');
        end

        % save vars
        if(mod(iter,var_save_invl)==0 || iter == 0)
            try
                save(sprintf('./postprocessing/phi_plot_%2d',iter),'phi_plot');
%                 save(sprintf('./postprocessing/theta_plot_%2d',iter),'theta_plot');
                save(sprintf('./postprocessing/tempr_plot_%2d',iter),'tempr_plot');
%                 save(sprintf('./postprocessing/conct_plot_%2d',iter),'conc_t_plot');
            catch
                fprintf('data write error skipped.\n');
            end
        end

    end
end

disp('All simulations complete!\n');