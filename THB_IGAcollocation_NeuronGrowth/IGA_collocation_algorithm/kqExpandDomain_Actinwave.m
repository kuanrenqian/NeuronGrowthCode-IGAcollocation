function [phi_out,theta_out,theta_ori_out,conc_t_out,param,tempr_out,phi_initial_out,theta_initial_out,tempr_initial_out,bcid] = ...
    kqExpandDomain_Actinwave(sz,phi,theta,conc_t,max_x,max_y,param,tempr,oldNuNv,newNuNv)

len = length(phi);
M = sqrt(len);

phi = reshape(phi,M,M);
theta = reshape(oldNuNv*theta,M,M);
conc_t = reshape(oldNuNv*conc_t,M,M);
actin = param.A;
hem = param.H;
HState = param.HState;
tempr = reshape(tempr,M,M);

% not necessary, just to be safe (stiffness setup insurance)
bcs = 4;
phi = phi(bcs:end-bcs,bcs:end-bcs);
theta = theta(bcs:end-bcs,bcs:end-bcs);
conc_t = conc_t(bcs:end-bcs,bcs:end-bcs);
actin = actin(bcs:end-bcs,bcs:end-bcs);
hem = hem(bcs:end-bcs,bcs:end-bcs);
HState = HState(bcs:end-bcs,bcs:end-bcs);
tempr = tempr(bcs:end-bcs,bcs:end-bcs);

[M, ~] = size(phi);
outSz = sqrt(sz);

phi_out =zeros(outSz);
theta_out =rand(outSz);
conc_t_out = zeros(outSz);
actin_out = zeros(outSz);
hem_out = zeros(outSz);
HState_out = zeros(outSz);
tempr_out =zeros(outSz);
% rot_map_out = zeros(outSz);

i_off = floor(outSz/2-M/2);
j_off = floor(outSz/2-M/2);

phi_out(3+i_off:M-2+i_off,3+j_off:M-2+j_off) = phi(3:M-2,3:M-2);
theta_out(3+i_off:M-2+i_off,3+j_off:M-2+j_off) = theta(3:M-2,3:M-2);
conc_t_out(3+i_off:M-2+i_off,3+j_off:M-2+j_off) = conc_t(3:M-2,3:M-2);

%%
theta_ori_out = theta_rotate_guide_sector(outSz,outSz,floor(M/2),floor(M/2),0);
%%

actin_out(3+i_off:M-2+i_off,3+j_off:M-2+j_off) = actin(3:M-2,3:M-2);
hem_out(3+i_off:M-2+i_off,3+j_off:M-2+j_off) = hem(3:M-2,3:M-2);
HState_out(3+i_off:M-2+i_off,3+j_off:M-2+j_off) = HState(3:M-2,3:M-2);
tempr_out(3+i_off:M-2+i_off,3+j_off:M-2+j_off) = tempr(3:M-2,3:M-2);
% rot_map_out(3+i_off:M-2+i_off,3+j_off:M-2+j_off) = rot_map(3:M-2,3:M-2);


phi_initial_out = reshape(phi_out,outSz,outSz);
theta_initial_out = reshape(theta_out,outSz,outSz);
tempr_initial_out = reshape(tempr_out,outSz,outSz);
for i = 2:outSz-1
    for j = 2:outSz-1
        phi_initial_out(i,j) = 0;
        theta_initial_out(i,j) = 0;
        tempr_initial_out(i,j) = 0;
    end
end
phi_initial_out = reshape(phi_initial_out,outSz*outSz,1);
theta_initial_out  = reshape(theta_initial_out,outSz*outSz,1);
tempr_initial_out  = reshape(tempr_initial_out,outSz*outSz,1);

phi_out =sparse(reshape(phi_out,outSz^2,1));
theta_out =sparse(newNuNv\reshape(theta_out,outSz^2,1));
conc_t_out =sparse(newNuNv\reshape(conc_t_out,outSz^2,1));
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

param.A = actin_out;
param.H = hem_out;
param.HState = HState_out;
