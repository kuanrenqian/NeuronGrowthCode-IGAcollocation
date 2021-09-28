function [colmats] = kqTHBders(Pixel,Jm,Pm,Dm)

NuNv = [];
N1uNv = [];
NuN1v = [];
N1uN1v = [];
N2uNv = [];
NuN2v = [];
lap = [];

[M,~] = size(Pixel);
Nx = sqrt(M);
Ny = Nx;

px = 0;
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

        for k = 1:ss
            CEb = Pm{SB(k,2),1};
            BB = Dm{SB(k,2),1};
            BB_active = BB{SB(k,1),10};
            
            NuNv(px,BB_active) = supp(k,1);
            N1uNv(px,BB_active) = suppx(k,1);
            NuN1v(px,BB_active) = suppy(k,1);
            N1uN1v(px,BB_active) = suppx(k,1) + suppy(k,1);
            N2uNv(px,BB_active) = suppxx(k,1);
            NuN2v(px,BB_active) = suppyy(k,1);
            lap(px,BB_active) = suppxx(k,1) + suppyy(k,1);
%             NuNv(BB_active,px) = supp(k,1);
%             N1uNv(BB_active,px) = suppx(k,1);
%             NuN1v(BB_active,px) = suppy(k,1);
%             N1uN1v(BB_active,px) = suppx(k,1) + suppy(k,1);
%             N2uNv(BB_active,px) = suppx(k,1);
%             NuN2v(BB_active,px) = suppyy(k,1);
%             lap(BB_active,px) = suppx(k,1) + suppyy(k,1);
        end
    end
end

colmats = struct('NuNv',NuNv,'N1uNv',N1uNv,'NuN1v',NuN1v,'N1uN1v', N1uN1v, ... 
    'N2uNv',N2uNv,'NuN2v',NuN2v,'lap',lap);
    