function [colmats] = kqTHBders(Pixel,Jm,Pm,Dm)

NuNv = [];
N1uNv = [];
NuN1v = [];
N1uN1v = [];
N2uNv = [];
NuN2v = [];
lap = [];

ac_ind_log = [];
BB_act_log = [];

[M,~] = size(Pixel);
Nx = sqrt(M);
Ny = Nx;

px = 0;
for i = 1:Nx
    for j = 1:Ny
        px = px +1;
        ac_ind = Pixel{px,1};
        ac_ind_log(end+1) = ac_ind;

        supp = Pixel{px,2};
        suppx = Pixel{px,3};
        suppy = Pixel{px,4};
        suppxx = Pixel{px,5};
        suppyy = Pixel{px,6};

        SB = Jm{ac_ind,1};
        ss = size(SB,1);

        for k = 1:ss
            CEb = Pm{SB(k,2),1};
            BB = Dm{SB(k,2),1};
            BB_active = BB{SB(k,1),10};
            BB_act_log(end+1) = BB_active;

            NuNv(px,BB_active) = supp(k,1);
            N1uNv(px,BB_active) = suppx(k,1);
            NuN1v(px,BB_active) = suppy(k,1);
            N1uN1v(px,BB_active) = suppx(k,1) + suppy(k,1);
            N2uNv(px,BB_active) = suppxx(k,1);
            NuN2v(px,BB_active) = suppyy(k,1);
            lap(px,BB_active) = suppxx(k,1) + suppyy(k,1);

%             %NuNv(NxNy,BB_active) = supp(k,1);
%             NuNv(ac_ind,BB_active) = supp(k,1);
%             N1uNv(ac_ind,BB_active) = suppx(k,1);
%             NuN1v(ac_ind,BB_active) = suppy(k,1);
%             N1uN1v(ac_ind,BB_active) = suppx(k,1) + suppy(k,1);
%             N2uNv(ac_ind,BB_active) = suppxx(k,1);
%             NuN2v(ac_ind,BB_active) = suppyy(k,1);
%             lap(ac_ind,BB_active) = suppxx(k,1) + suppyy(k,1);

        end
    end
end

colmats = struct('NuNv',NuNv,'N1uNv',N1uNv,'NuN1v',NuN1v,'N1uN1v', N1uN1v, ... 
    'N2uNv',N2uNv,'NuN2v',NuN2v,'lap',lap,'BB_act_log',BB_act_log);
    