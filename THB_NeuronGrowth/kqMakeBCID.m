function [bcid,bcid_cp] = kqMakeBCID(bc_sz,bc_layers,coll_X,coll_Y,THBfinal)

bcid = zeros([bc_sz,bc_sz]);
for i = 1:bc_sz
    bcid(1:1+bc_layers,i) = ones(bc_layers+1,1);
    bcid(bc_sz-bc_layers:bc_sz,i) = ones(bc_layers+1,1);
    bcid(i,1:1+bc_layers) = ones(1,bc_layers+1);
    bcid(i,bc_sz-bc_layers:bc_sz) = ones(1,bc_layers+1);
end

bcid_cp  = interp2(coll_X,coll_Y,bcid,THBfinal(:,1),THBfinal(:,2),'spline');
bcid_cp(isnan(bcid_cp(:))) = 0.0;
bcid_cp = round(bcid_cp);

% bcid_cp = zeros(length(THBfinal),1);
% for i = 1:length(THBfinal)
%     x = THBfinal(i,1);
%     y = THBfinal(i,2);
%     if (x==0 || y ==0 || x==1 || y==1)
%         bcid_cp(i) = 1;
%     end
% end

end