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

M = length(Pixel);

for i = 1:M
    ac_ind = Pixel{i,1};
    ac_ind_log(end+1) = ac_ind;
    
    supp = Pixel{i,2};
    suppx = Pixel{i,3};
    suppy = Pixel{i,4};
    suppxx = Pixel{i,5};
    suppyy = Pixel{i,6};

    SB = Jm{ac_ind,1};
    ss = size(SB,1);

    for k = 1:ss
        CEb = Pm{SB(k,2),1};
        BB = Dm{SB(k,2),1};
        BB_active = BB{SB(k,1),10};
        BB_act_log(end+1) = BB_active;

        NuNv(i,BB_active) = supp(k,1);
        N1uNv(i,BB_active) = suppx(k,1);
        NuN1v(i,BB_active) = suppy(k,1);
        N1uN1v(i,BB_active) = suppx(k,1) + suppy(k,1);
        N2uNv(i,BB_active) = suppxx(k,1);
        NuN2v(i,BB_active) = suppyy(k,1);
        lap(i,BB_active) = suppxx(k,1) + suppyy(k,1);

    end
end
    
% store the LU decomposition of NuNv
[L_NuNv, U_NuNv] = lu(NuNv);

% store the diagonals of NuNv, N1uNv, NuN1v
[NuNv_flip, NuNv_id] = extract_diags(NuNv);
[N1uNv_flip, N1uNv_id] = extract_diags(N1uNv);
[NuN1v_flip, NuN1v_id] = extract_diags(NuN1v);
colmats = struct('NuNv', NuNv, 'NuNv_flip', NuNv_flip, 'NuNv_id', NuNv_id,...
        'L_NuNv',L_NuNv,'U_NuNv',U_NuNv,'N1uNv',N1uNv,'N1uNv_flip',N1uNv_flip,...
        'N1uNv_id',N1uNv_id,'NuN1v',N1uNv,'NuN1v_flip',NuN1v_flip,'NuN1v_id',...
        NuN1v_id,'N1uN1v', N1uN1v, 'N2uNv',N2uNv,'NuN2v',NuN2v,'lap',lap);

% colmats = struct('NuNv',NuNv,'N1uNv',N1uNv,'NuN1v',NuN1v,'N1uN1v', N1uN1v, ... 
%     'N2uNv',N2uNv,'NuN2v',NuN2v,'lap',lap,'ac_ind_log',ac_ind_log,'BB_act_log',BB_act_log);