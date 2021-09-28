function [out_cp] = kqConvert2refinedPts(input,CEb,THBfinal)

[M,~] = size(CEb);
coll_sz = sqrt(M);

% create collocation mesh grid
coll_X_array = CEb(1:coll_sz,2).';
coll_Y_array = CEb(1:coll_sz,2).';
[coll_X,coll_Y] = meshgrid(coll_X_array,coll_Y_array);

input = reshape(input,coll_sz,coll_sz);
out_cp  = interp2(coll_X,coll_Y,input,THBfinal(:,1),THBfinal(:,2),'spline');
out_cp(isnan(out_cp(:))) = 0.0;

end