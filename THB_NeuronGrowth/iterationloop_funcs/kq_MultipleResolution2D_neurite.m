iteration=0;
iter_vector = [];
n_vector = [];
pU = parameters.pU;
pV = parameters.pV;

%% Multilevel framework
for level = parameters.maxlevel:-1:1
        
    disp(['Registration level:' num2str(parameters.maxlevel-level+1)]);
    
    %% Meshgrid for later use (will be used to change size in the future)
    nx = Nx;
    ny = Ny;
    [X,Y] = meshgrid(linspace(0,1,Nx),linspace(0,1,Ny));

    %% Construct B-spline grid
    maxlev = parameters.maxlevel-level+1;
    [Dm,Pm,Em,Bvect,knotvectorU,knotvectorV,nobU,nobV,nelemU] = setBsplineGrid(maxlev,parameters,Nx,Ny,dx,dy);
    
    for multilev = 0:1:maxlev-1
        if(multilev>0)
            % local refinement based on laplacian of phi
            for j =1:floor(bf_ct/2)
%                 if(rf_cp(j)~=0) % laplacian of phi
                    bbc = bf(j,1:2);
                    bf_lev = bf(j,3);
%                     [Dm,Em,Pm] =  Refine2D(bbc(1,1),bbc(1,2),bf_lev,Dm,Em,Pm,knotvectorU,knotvectorV,pU,pV);
                    [Dm,Em,Pm] =  Refine2Dtrunc1(bbc(1,1),bbc(1,2),bf_lev,Dm,Em,Pm,knotvectorU,knotvectorV,pU,pV);
%                 end
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
    end
        
    % rf_cp is purely for refinement, actual phi is initialized on locally
    % refined points
    Pfinal = zeros(bf_ct,2);
    for i = 1:bf_ct
        bbc = bf(i,1:2);
        bf_lev = bf(i,3);
        bi = nobU(bf_lev,1)*(bbc(1,2)-1)+bbc(1,1);
        Pi = Pm{bf_lev,1};
        Pfinal(i,1) = Pi(bi,1);
        Pfinal(i,2) = Pi(bi,2);
    end
    [LAP] = del2(reshape(phi,Nx,Ny));
    rf_cp  = interp2(X,Y,LAP,Pfinal(:,2),Pfinal(:,1));
    rf_cp(isnan(rf_cp(:))) = 0;

    %if level == 1 % second level, 1 local refinement
        figure;
        displayAdaptiveGrid(ac,Coeff,Em,knotvectorU,knotvectorV,Jm,Pm,parameters,dx*nx,dy*ny);
        title(sprintf('Mesh with %d refinements',maxlev-1));
        axis square;
        drawnow;
        kq_iterationloop_neurite
    %end
end