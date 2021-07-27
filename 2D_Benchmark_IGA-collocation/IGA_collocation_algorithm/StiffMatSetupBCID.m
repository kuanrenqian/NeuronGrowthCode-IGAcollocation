function [coll_Lhs, coll_Rhs] = StiffMatSetupBCID(coll_Lhs, coll_Rhs,bcid,N)

for i=1:length(bcid)
    % change stiffness mat where there is non-zero boundary condition
    if bcid(i)==1
        coll_Lhs(i,:) = zeros(1,length(bcid));
        coll_Lhs(:,i) = zeros(1,length(bcid))';
        coll_Lhs(i,i) = 1;
        coll_Rhs(i) = N(i);
    end
end
coll_Lhs = sparse(coll_Lhs);
coll_Rhs = sparse(coll_Rhs);