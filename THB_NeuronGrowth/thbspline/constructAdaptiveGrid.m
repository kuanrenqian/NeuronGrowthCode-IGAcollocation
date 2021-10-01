function [Jm,Coeff,Pixel] = constructAdaptiveGrid(ac,parameters,Dm,Em,X,Y,knotvectorU,knotvectorV,multilev,nobU,nobV,nelemU)

ac_ct = size(ac,1);
maxlev=size(knotvectorU,1);
Jm = cell(ac_ct,1);

pU = parameters.pU;
pV = parameters.pV;

for i = 1:ac_ct
    
    cell_ind = ac(i,1);
    cell_lev = ac(i,2);
    counter = 0;
    DE = Dm{cell_lev,1};
    EE = Em{cell_lev,1};
    
    %Collect the non-zeros HB splines over each active cell
    local_b = EE{cell_ind,6};
    lb_size = size(local_b,2);
    cell_supp = zeros(1,2);
    
    for j = 1:lb_size
        if(DE{local_b(1,j),3}==1)
            counter= counter+1;
            cell_supp(counter,1) = local_b(1,j);
            cell_supp(counter,2) = cell_lev;
        end
    end
    curr_cell = cell_ind;
    
    if(cell_lev~=1)
        for m = (cell_lev-1):-1:1
            EEM = Em{m,1};
            EEMM = Em{m+1,1};
            DEM = Dm{m,1};
            
            curr_cell = EEMM{curr_cell,10};
            local_supp = EEM{curr_cell,6};
            local_supp_size = size(local_supp,2);
            for j = 1:local_supp_size
                if(DEM{local_supp(1,j),3}==1)
                    counter=counter+1;
                    cell_supp(counter,1) = local_supp(1,j);
                    cell_supp(counter,2) = m;
                end
            end
        end
    end
    Jm{i,1} = cell_supp;
end
%==========================================================================
Coeff = cell(ac_ct,1);
for  i = 1:ac_ct
    count=0;
    cell_ind = ac(i,1);
    cell_lev = ac(i,2);
    ctemp = cell(maxlev,1);
    EE = Em{cell_lev,1};
    DE = Dm{cell_lev,1};
    local_b = EE{cell_ind,6};
    lb_size = size(local_b,2);
    ct = zeros(16,16);
    ct11 = zeros(16,16);
    coef_arr = zeros(1,16);
    for j = 1:lb_size
        ct(j,j) = 1;
    end
    ctemp{cell_lev,1} = ct;
    for j =1:lb_size
        if(DE{local_b(1,j),3}==1)
            count= count+1;
            coef_arr(count,:) = ct(j,:);
        end
    end
    curr_cell = cell_ind;
    curr_cell_level = cell_lev;
    if(cell_lev~=1)
        for m=(cell_lev-1):-1:1
            EEM = Em{m,1};
            DEM = Dm{m,1};
            DEMM = Dm{m+1,1};
            EEMM = Em{m+1,1};
            local_s1 = EEMM{curr_cell,6};
            curr_cell = EEMM{curr_cell,10};
            local_s = EEM{curr_cell,6};
            ls_size = size(local_s,2);
            ct = zeros(16,16);
            ct11 = zeros(16,16);
            for j=1:ls_size
                cmat = DEM{local_s(1,j),12};
                for k=1:ls_size
                    for k1=1:size(cmat,1)
                        if(cmat(k1,2)==local_s1(1,k))
                            ct(j,k) = cmat(k1,1);
                            if(DEMM{local_s1(1,k),3}==0&&DEMM{local_s1(1,k),7}==0)
                                ct11(j,k) = ct(j,k);
                            end
                        end
                    end
                end
            end
            ct1 = ctemp{m+1,1};
            ct = ct*ct1;
            ctnew = ct11*ct1;
            ctemp{m,1} = ct;
            for j=1:ls_size
                if(DEM{local_s(1,j),3}==1)
                    count = count+1;
                    coef_arr(count,:) = ctnew(j,:);
                end
            end
        end
    end
    Coeff{i,1} = coef_arr;
end

Pixel = cell(size(X,1)*size(X,2),2);
kuMAX = knotvectorU{multilev+1,1};
kvMAX = knotvectorV{multilev+1,1};
unobMAX = nobU(multilev+1,1);
vnobMAX = nobV(multilev+1,1);
EE = Em{multilev+1,1};
p_ind = 0;
for i = 1:size(X,1)
    for j = 1:size(X,2)
        p_ind = p_ind+1;
        uu = FindSpan(unobMAX-1,pU,X(j,i),kuMAX) + 1;
        vv = FindSpan(vnobMAX-1,pV,Y(j,i),kvMAX) + 1;
        
        cellx = (uu-pU);
        celly = (vv-pV);
        cell_ind = nelemU(multilev+1,1)*(celly-1)+cellx;
        curr_cell = cell_ind;
        if(EE{cell_ind,4} == 1)
            
            act_ind = EE{cell_ind,11};
            SB = Jm{act_ind,1};
            sb_size = size(SB,1);
            pix_coeff = Coeff{act_ind,1};
            SB1 = EE{cell_ind,6};
            sb_size1 = size(SB1,1);
            
            RRD1 = dersbasisfuns3(uu,pU,2,X(j,i),kuMAX);
            RRD2 = dersbasisfuns3(vv,pV,2,Y(j,i),kvMAX);

            RRD = (RRD1(1,:)')*(RRD2(1,:));
            RRDU = (RRD1(1,:)')*(RRD2(2,:));
            RRDV = (RRD1(2,:)')*(RRD2(1,:));
            RRDUU = (RRD1(1,:)')*(RRD2(3,:));
            RRDVV = (RRD1(3,:)')*(RRD2(1,:));
            
            inc=0;
            phii = zeros(16,1);
            phiiu = zeros(16,1);
            phiiv = zeros(16,1);
            phiiuu = zeros(16,1);
            phiivv = zeros(16,1);
            
            for m1 = 1:size(RRD,1)
                for m2=1:size(RRD,2)
                    phii(16-inc,1)= RRD(m2,m1);
                    phiiu(16-inc,1) = RRDU(m2,m1);
                    phiiv(16-inc,1) = RRDV(m2,m1);
                    phiiuu(16-inc,1) = RRDUU(m2,m1);
                    phiivv(16-inc,1) = RRDVV(m2,m1);
                    inc=inc+1;
                end
            end
            
            phi_pi = pix_coeff*phii;
            phi_piu = pix_coeff*phiiu;
            phi_piv = pix_coeff*phiiv;
            phi_piuu = pix_coeff*phiiuu;
            phi_pivv = pix_coeff*phiivv;
            
            Pixel{p_ind,1} = act_ind;
            Pixel{p_ind,2} = phi_pi;
            Pixel{p_ind,3} = phi_piu;
            Pixel{p_ind,4} = phi_piv;
            Pixel{p_ind,5} = phi_piuu;
            Pixel{p_ind,6} = phi_pivv;
        else
            
            for m = (multilev+1):-1:2
                EEM = Em{m,1};
                EEMM = Em{(m-1),1};
                curr_cell = EEM{curr_cell,10};
                if(EEMM{curr_cell,4} == 1)
                    
                    act_ind = EEMM{curr_cell,11};
                    SB = Jm{act_ind,1};
                    sb_size = size(SB,1);
                    knotuu = knotvectorU{(m-1),1};
                    knotvv = knotvectorV{(m-1),1};
                    unob1 = nobU(m-1,1);
                    vnob1 = nobV(m-1,1);
                    u1 = FindSpan(unob1-1,pU,X(j,i),knotuu) + 1;
                    v1 = FindSpan(vnob1-1,pV,Y(j,i),knotvv) + 1;
                    pix_coeff = Coeff{act_ind,1};
                    SB1 = EE{cell_ind,6};
                    sb_size1 = size(SB1,1);

                    RRD1 = dersbasisfuns3(u1,pU,2,X(j,i),knotuu);
                    RRD2 = dersbasisfuns3(v1,pV,2,Y(j,i),knotvv);

                    RRD = (RRD1(1,:)')*(RRD2(1,:));
                    RRDU = (RRD1(1,:)')*(RRD2(2,:));
                    RRDV = (RRD1(2,:)')*(RRD2(1,:));
                    RRDUU = (RRD1(3,:)')*(RRD2(1,:));
                    RRDVV = (RRD1(1,:)')*(RRD2(3,:));
                    
                    inc=0;
                    phii = zeros(16,1);
                    phiiu = zeros(16,1);
                    phiiv = zeros(16,1);
                    phiiuu = zeros(16,1);
                    phiivv = zeros(16,1);
                    for m1 = 1:size(RRD,1)
                        for m2=1:size(RRD,2)
                            phii(16-inc,1)= RRD(m2,m1);
                            phiiu(16-inc,1) = RRDU(m2,m1);
                            phiiv(16-inc,1) = RRDV(m2,m1);
                            phiiuu(16-inc,1) = RRDUU(m2,m1);
                            phiivv(16-inc,1) = RRDVV(m2,m1);
                            inc=inc+1;
                        end
                    end
                    phi_pi = pix_coeff*phii;
                    phi_piu = pix_coeff*phiiu;
                    phi_piv = pix_coeff*phiiv;
                    phi_piuu = pix_coeff*phiiuu;
                    phi_pivv = pix_coeff*phiivv;
                    
                    Pixel{p_ind,1} = act_ind;
                    Pixel{p_ind,2} = phi_pi;
                    Pixel{p_ind,3} = phi_piu;
                    Pixel{p_ind,4} = phi_piv;
                    Pixel{p_ind,5} = phi_piuu;
                    Pixel{p_ind,6} = phi_pivv;
                    break;
                end
            end
        end
    end
end

end

