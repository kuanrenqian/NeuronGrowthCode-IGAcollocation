function param = GetParam(phi,dt)
% load('phiplot_05202021');
% load('phiplot_05202021_100by100');
len = length(phi);
phi = full(phi);
p = phi;
for i=1:len^2
    if(phi(i)>0.5)
        p(i) = 1;
    else
        p(i) = 0;
    end
end
param.RegionMask = p;
param.phi = reshape(phi,len^2,1);

% 0) Make simulations deterministic
param.state = 136;
rng(param.state);

% 1) Domain parameters
% T:  total simulation time T 
% W:  width W of disk domain
% dt: H nucleation updated every time step dt; 
%     A, H, R updated continuously in [t,t+dt]
param.T  = 1000;
param.dt = dt;
param.T_current = 0;
param.W = floor((len-3)/2)+1;
param.N = 2*param.W +1;

% 2) Model parameters
%  (a) Kernel [G_SA] = Short range Activator
%   sigma: steepness of Gaussian, radius: size of kernel window. 
%   radius = 0: Gaussian assumed so wide that it's constant = 1/(area disk)
param.kernel.sigma  = 1;
param.kernel.radius = 1; 
param.kernel.width  = 2 * param.kernel.radius +1;
param.G_SA = fspecial('gaussian',param.kernel.width(1),param.kernel.sigma(1));

param.c_A  = 1;
param.c_H  = [1 0.05];
param.HTotal = 2000;
param.H_thresh_on  = 0.1;
param.H_thresh_off = 0.01;

% 3) Initialize structures 
param.A_0  = zeros(param.N);
param.H_0  = zeros(param.N);
param.A    = zeros(param.N);
param.H    = zeros(param.N);

% state 0=off 1=on
param.HState = zeros(param.N);

% 4) Set up initial conditions
param.A_0(param.W-1:param.W+1,param.W-1:param.W+1) = 0.1;
param.H_0(param.W-1:param.W+1,param.W-1:param.W+1) = 0.5;
param.HState(param.W-1:param.W+3,param.W-1:param.W+3) = 1;

param.version = 'Hem1_v11.m';
param.A = param.A_0;
param.H = param.H_0;
end
