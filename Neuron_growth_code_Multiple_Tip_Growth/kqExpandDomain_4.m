function [phi_out,theta_out,conc_out,tempr_out,phi_initial_out,theta_initial_out,conc_initial_out,tempr_initial_out,bcid] = ...
    kqExpandDomain_4(sz,phi,theta,conc,tempr,oldNuNv,newNuNv)

len = length(phi);
M = sqrt(len);
phi = reshape(phi,M,M);
theta = reshape(oldNuNv*theta,M,M);
conc = reshape(conc,M,M);
tempr = reshape(tempr,M,M);

phi = phi(5:end-5,5:end-5);
theta = theta(5:end-5,5:end-5);
conc = conc(5:end-5,5:end-5);
tempr = tempr(5:end-5,5:end-5);

[M, ~] = size(phi);
outSz = sqrt(sz);

phi_out =zeros(outSz);
theta_out =rand(outSz);
for i = 1:outSz
    for j = 1:outSz
%         theta_out(i,j) = i/outSz*pi;

        x=i-outSz/2;
        y=j-outSz/2;
        theta_out(i,j) = ((atan2(y,x)));
%         theta_out(i,j) = round(2*(atan2(y,x)))/2;
        if isnan(theta_out(i,j))
            theta_out(i,j) = 1;
        end
%         theta_out(i,j) = theta_out(i,j);
%         if(theta_out(i,j)>pi)
%             theta_out(i,j) = theta_out(i,j) - 2*pi;
%         elseif (theta_out(i,j)<-pi)
%             theta_out(i,j) = theta_out(i,j) + 2*pi;
%         end
        
    end
end

%%
% for i = 1:length(theta_out)
% %     if (phi(i)<1&&j>max_y)
%         theta_out(i) = theta_out(i)+rotate;
%         if(theta_out(i)>pi)
%             theta_out(i) = theta_out(i) - 2*pi;
%         elseif (theta_out(i)<-pi)
%             theta_out(i) = theta_out(i) + 2*pi;
%         end
% %     end
% end

    %%
conc_out =ones(outSz).*0.75;
tempr_out =zeros(outSz);

i_off = floor(outSz/2-M/2);
j_off = floor(outSz/2-M/2);

for i = 3:M-2
    for j = 3:M-2
        phi_out(i+i_off, j+j_off) = phi(i,j);
    end
end
for i = 3:M-2
    for j = 3:M-2
        theta_out(i+i_off, j+j_off) = theta(i,j);
    end
end
for i = 3:M-2
    for j = 3:M-2
        conc_out(i+i_off, j+j_off) = conc(i,j);
    end
end
for i = 3:M-2
    for j = 3:M-2
        tempr_out(i+i_off, j+j_off) = tempr(i,j);
    end
end

phi_initial_out = reshape(phi_out,outSz,outSz);
conc_initial_out = reshape(conc_out,outSz,outSz);
theta_initial_out = reshape(theta_out,outSz,outSz);
tempr_initial_out = reshape(tempr_out,outSz,outSz);
for i = 2:outSz-1
    for j = 2:outSz-1
        phi_initial_out(i,j) = 0;
        conc_initial_out(i,j) = 0;
        theta_initial_out(i,j) = 0;
        tempr_initial_out(i,j) = 0;
    end
end
phi_initial_out = reshape(phi_initial_out,outSz*outSz,1);
conc_initial_out = reshape(conc_initial_out,outSz*outSz,1);
theta_initial_out  = reshape(theta_initial_out,outSz*outSz,1);
tempr_initial_out  = reshape(tempr_initial_out,outSz*outSz,1);

phi_out =sparse(reshape(phi_out,outSz^2,1));
theta_out =sparse(newNuNv\reshape(theta_out,outSz^2,1));
conc_out =sparse(newNuNv\reshape(conc_out,outSz^2,1));
tempr_out =sparse(reshape(tempr_out,outSz^2,1));

bcid = zeros([outSz,outSz]);
for i = 1:outSz
    bcid(1,i) = 1;
    bcid(end,i) = 1;
    bcid(i,1) = 1;
    bcid(i,end) = 1;
end
bcid = reshape(bcid,outSz^2,1);
bcid = sparse(bcid);
