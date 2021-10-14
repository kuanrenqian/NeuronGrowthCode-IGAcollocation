function [bcid,bcid_cp] = kqMakeBCID(bc_sz,bc_layers,X,Y,THBfinal)

% build boundary condition id as a binary field
bcid = zeros([bc_sz,bc_sz]);
for i = 1:bc_sz
    bcid(1:bc_layers,i) = ones(bc_layers,1);
    bcid(bc_sz-bc_layers+1:bc_sz,i) = ones(bc_layers,1);
    bcid(i,1:bc_layers) = ones(1,bc_layers);
    bcid(i,bc_sz-bc_layers+1:bc_sz) = ones(1,bc_layers);
end

% interp to THBfinal points
bcid_cp  = interp2(X,Y,bcid,THBfinal(:,1),THBfinal(:,2),'spline');
bcid_cp(isnan(bcid_cp(:))) = 0.0;
bcid_cp = round(bcid_cp);
end