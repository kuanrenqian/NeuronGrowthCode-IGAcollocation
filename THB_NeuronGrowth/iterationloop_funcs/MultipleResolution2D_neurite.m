%%
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
        
        [X,Y] = meshgrid(linspace(0,1,Nx*2^(parameters.maxlevel-level)),linspace(0,1,Ny*2^(parameters.maxlevel-level)));
    else
        phi_new = zeros((Nx)* 2  ,(Ny)* 2  );
        
        phi_new(floor((Nx* 2  -Nx)/2)+1:Nx+floor((Nx* 2  -Nx)/2),floor((Ny* 2  -Ny)/2)+1:Ny+floor((Ny* 2  -Ny)/2)) = phi;

        Nx = Nx*2;
        Ny = Ny*2;
        NxNy = (Nx)*(Ny);
                
        phi = phi_new;
        
        nx = Nx;
        ny = Ny;

       [X,Y] = meshgrid(linspace(0,1,Nx),linspace(0,1,Ny));
    end
    
    
    %% Construct B-spline grid
    maxlev = parameters.maxlevel-level+1;
    [Dm,Pm,Em,Bvect,knotvectorU,knotvectorV,nobU,nobV,nelemU] = setBsplineGrid(maxlev,parameters,Nx,Ny,dx,dy);
    
    for multilev = 0:1:maxlev-1
        if(multilev>0)
%             grad_log = [];
            for j =1:floor(bf_ct/2)
                bbc = bf(j,1:2);
                bf_lev = bf(j,3);
                [Dm,Em,Pm] =  Refine2Dtrunc1(bbc(1,1),bbc(1,2),bf_lev,Dm,Em,Pm,knotvectorU,knotvectorV,pU,pV);
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
        
        phi_cp  = interp2(X,Y,reshape(phi,Nx,Ny),Pfinal(:,2),Pfinal(:,1),'spline');
        phi_cp(isnan(phi_cp(:))) = 0.0;
        cell_co = zeros(ac_ct,2);
        for i = 1:ac_ct
            cell_id = ac(i,1);
            cell_le = ac(i,2);
            EEM = Em{cell_le,1};
            cell_co(i,1) = EEM{cell_id,8};
            cell_co(i,2) = EEM{cell_id,9};
        end
        
        phi_imgg = reshape(phi,Nx,Ny).*255;
        [phidxi, phidyi] = gradient(phi_imgg);
        Idiff = sqrt(phidxi.^2+phidyi.^2);
        Cell_grad = interp2(X,Y,Idiff,cell_co(:,2),cell_co(:,1),'spline');
        meangrad = mean2(Idiff);
    end

    kq_laplace_test
end

    