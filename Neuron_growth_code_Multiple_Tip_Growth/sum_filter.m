function [phi_sum] = sum_filter(phi,state)
% THis function takes phi as input and output a phi_sum variable that has
% maximum value at tips (0~1)

% get size of input
[Nx,Ny] = size(phi);
% round phi -> discrete
phi = round(phi);
% initialize
phi_sum = zeros(Nx,Ny);
% loop through and calculate sum of phi values around i,j
for i = 5:Nx-4
    for k = 5:Ny-4
        for j = k-4:k+4
            phi_sum(i,k) = phi_sum(i,k) + sum(phi(i-4,j) + phi(i-3,j) + ...
                phi(i-2,j) + phi(i-1,j) + phi(i,j) + phi(i+1,j) + ...
                phi(i+2,j) + phi(i+3,j) + phi(i+4,j));
        end
    end
end

%% remember to simplify (using imlocalmin, then this section is not needed)
% scaling for better identification
phi_sum = phi./phi_sum;
phi_sum_max = max(max(phi_sum));
phi_sum = phi_sum./phi_sum_max;
phi_sum(isnan(phi_sum))=0;

phi_sum(phi_sum<0.7)=0;

% if state == 1
%     for i = 1:Nx
%         for j = 1:Ny
%             r = sqrt((i-Nx/2)^2+(j-Ny/2)^2);
%             if r<18
%                 phi_sum(i,j) = 0;
%             end
%         end
%     end
% end