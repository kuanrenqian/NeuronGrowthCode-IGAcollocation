function [lhs, rhs] = kqNGStiffMatSetupBCID(lhs_initial,rhs_initial,initial_cp,bcid_cp)

lhs = lhs_initial(:,bcid_cp~=1);
C = lhs_initial(:,bcid_cp==1);
% requires setting up initial_cp correctly outside this fun
f = initial_cp(bcid_cp==1); rhs = rhs_initial-C*f;

%% Backup
% [row,col] = size(coll_Lhs);
% indL = [];
% for i=1:col
%     if bcid_cp(i)==1
%         coll_Lhs(:,i) = zeros(1,row)';
%     end
% end
