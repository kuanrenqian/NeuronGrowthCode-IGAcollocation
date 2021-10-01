function [coll_Lhs, coll_Rhs] = kqStiffMatSetupBCID(coll_Lhs,coll_Rhs,bcid)

[row,col] = size(coll_Lhs);
% ind = [];
for i=1:col
    if bcid(i)==1
        coll_Lhs(:,i) = zeros(1,row)';
%         ind(end+1) = i;
    end
end

% coll_Lhs(:,ind) = [];