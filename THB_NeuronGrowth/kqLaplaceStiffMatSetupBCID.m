function [lhs, rhs] = kqLaplaceStiffMatSetupBCID(lap,initial_cp,bcid_cp)

lhs = lap(:,bcid_cp~=1);
C = lap(:,bcid_cp==1);

% requires setting up initial_cp correctly outside this fun
f = initial_cp(bcid_cp==1); 
rhs = -C*f;

%% Backup
% [row,col] = size(coll_Lhs);
% indL = [];
% for i=1:col
%     if bcid_cp(i)==1
%         coll_Lhs(:,i) = zeros(1,row)';
%     end
% end
