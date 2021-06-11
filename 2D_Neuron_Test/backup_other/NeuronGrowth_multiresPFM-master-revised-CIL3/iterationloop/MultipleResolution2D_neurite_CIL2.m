%function [phi,tempr] = MultipleResolution2D_pf(M,F,phi,parameters)

iteration=0;
iter_vector = [];
n_vector = [];
pU = parameters.pU;
pV = parameters.pV;

%% Multilevel framework
for level = parameters.maxlevel:-1:1
    
    disp(['Registration level:' num2str(parameters.maxlevel-level+1)]);
    
    %% downsample image
    if(level==parameters.maxlevel)
        
        nx = Nx;
        ny = Ny;
        
        [X,Y] = meshgrid(dx:dx:dx*(Nx),dy:dy:(Ny)*dy);
        [laplacian_mat] = laplacian(Nx,Ny,dx,dy);
        flag_energy = zeros(Nx,Ny);
    else
        phi_new = zeros((Nx)* 2  ,(Ny)* 2  );
        phi_new3 = zeros((Nx)* 2  ,(Ny)* 2  );
        phi_ori_new = zeros((Nx)* 2  ,(Ny)* 2  );
        tempr_new = zeros((Nx)* 2  ,(Ny)* 2  );
        conc_new = 0.75.*ones((Nx)* 2,(Ny)* 2  );
        theta_new = rand((Nx)* 2  ,(Ny)* 2  );
        theta_new_original = theta_new;
        conct_new = zeros((Nx)*2,(Ny)*2);
        flag_energy_new = zeros((Nx)*2,(Ny)*2);
        
        phi_new(floor((Nx* 2  -Nx)/2)+1:Nx+floor((Nx* 2  -Nx)/2),floor((Ny* 2  -Ny)/2)+1:Ny+floor((Ny* 2  -Ny)/2)) = phi;
        phi_new3(floor((Nx* 2  -Nx)/2)+1:Nx+floor((Nx* 2  -Nx)/2),floor((Ny* 2  -Ny)/2)+1:Ny+floor((Ny* 2  -Ny)/2)) = phi3;
        phi_ori_new(floor((Nx* 2  -Nx)/2)+1:Nx+floor((Nx* 2  -Nx)/2),floor((Ny* 2  -Ny)/2)+1:Ny+floor((Ny* 2  -Ny)/2)) = phi_original;
        
        
        tempr_new(floor((Nx* 2  -Nx)/2)+1:Nx+floor((Nx* 2  -Nx)/2),floor((Ny* 2  -Ny)/2)+1:Ny+floor((Ny* 2  -Ny)/2)) = tempr;
        theta_new(floor((Nx* 2  -Nx)/2)+1:Nx+floor((Nx* 2  -Nx)/2),floor((Ny* 2  -Ny)/2)+1:Ny+floor((Ny* 2  -Ny)/2)) = theta;
        theta_new_original(floor((Nx* 2  -Nx)/2)+1:Nx+floor((Nx* 2  -Nx)/2),floor((Ny* 2  -Ny)/2)+1:Ny+floor((Ny* 2  -Ny)/2)) = theta_original;
        conc_new(floor((Nx* 2  -Nx)/2)+1:Nx+floor((Nx* 2  -Nx)/2),floor((Ny* 2  -Ny)/2)+1:Ny+floor((Ny* 2  -Ny)/2)) = conc;
        conct_new(floor((Nx*2-Nx)/2)+1:Nx+floor((Nx*2-Nx)/2),floor((Ny*2-Ny)/2)+1:Ny+floor((Ny*2-Ny)/2)) = conc_t;
        flag_energy_new(floor((Nx*2-Nx)/2)+1:Nx+floor((Nx*2-Nx)/2),floor((Ny*2-Ny)/2)+1:Ny+floor((Ny*2-Ny)/2)) = flag_energy;
        
        Nx = Nx*2;
        Ny = Ny*2;
        NxNy = (Nx)*(Ny);
        
        [laplacian_mat] = laplacian(Nx,Ny,dx,dy);
        
        phi = phi_new;
        phi3 = phi_new3;
        phi_original = phi_ori_new;
        tempr = tempr_new;
        theta = theta_new;
        theta_original = theta_new_original;
        conc = conc_new;
        conc_t = conct_new;
        flag_energy = flag_energy_new;
        flag_center = zeros(Nx,Ny);
        for i=1:Nx
            for j=1:Ny
                if(j<=Nx/2)
                    flag_center(i,j) = 1;
                end
            end
        end
        
        nx = Nx;
        ny = Ny;
        
        [X,Y] = meshgrid(dx:dx:(dx*Nx),dy:dy:(Ny)*dy);
        
        if(level==2)
            nstep = 1900;
        end
        if(level==1)
            nstep=2000;
        end
        
        center(:,1) = center(:,1)+(Nx-Nx/2)/2;
        center(:,2) = center(:,2)+(Ny-Ny/2)/2;
    end
    
    %% Pre-allocating variables
    px = 0;
    phi1 = zeros(Nx,Ny);
    phidx1= zeros(Nx,Ny);
    phidy1= zeros(Nx,Ny);
    
    NxNy = Nx*Ny;
    ac_ind = zeros(NxNy,1);
    suppx = cell(NxNy,1);
    suppy = cell(NxNy,1);
    Ind = cell(NxNy,1);
    FFX = zeros(1,NxNy);
    FFY = zeros(1,NxNy);
    
    %% Construct B-spline grid
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
    
    if(level==1)
        displayAdaptiveGrid(ac,Coeff,Em,knotvectorU,knotvectorV,Jm,Pm,parameters,dx*nx,dy*ny);
    end
    
    [phidyo,phidxo] = gradient_mat(phi_original,Nx,Ny,dx,dy);
    sq_grad_phi = phidyo.^2+phidxo.^2;
    sum_grad_phi_ori = sum(sq_grad_phi(:));
    s_coeff = 3.e-6.*ones(Nx,Ny);
    
    save('theta_CIL2_original_new.mat','theta_original');
    
    iterationloop_neurite_CIL2
end