close all;
clear;
clc;

addpath('../IGA_collocation_algorithm');

% 1) get simulation parameters
param = GetParam;

%% Start video writer
% writer = VideoWriter('active_wave_05262021'); % Define VideoWriter object
% open(writer);

%%
figure(1);
set(gcf, 'position',[100 100 1200 400]);

for t = 1:param.T
    disp(t);

    % 2) set up information for ODE solver
    T   = [0 param.dt];
    sz  = size(param.A);
    
%     collA = reshape(NuNv*reshape(param.A,sz(1)*sz(2),1),sz);
%     collH = reshape(NuNv*reshape(param.H,sz(1)*sz(2),1),sz);
    
    Y_0 = [reshape(param.A,sz(1)*sz(2),1); reshape(param.H,sz(1)*sz(2),1)];
%     Y_0 = [reshape(collA,sz(1)*sz(2),1); reshape(collH,sz(1)*sz(2),1)];
    
    % 3) call ODE solver
    [~,Y] = ode45(@evolvesystem,T,Y_0,[],param);

    sz = size(param.H);
    param.A = reshape(Y(end,1:end/2),sz);
    param.H = reshape(Y(end,end/2+1:end),sz);
%     collA = reshape(Y(end,1:end/2),sz);
%     collH = reshape(Y(end,end/2+1:end),sz);
    param.T_current = t*param.dt;

    %% 4) update Hem
    % 1) Change from state 0 to 1
    H_C = max(0,(param.HTotal - sum(sum(param.H)))/param.HTotal);
    kSA = conv2(param.H,param.G_SA,'same');
    lambda = H_C*param.dt*kSA;

    Q = poissrnd(lambda,size(param.H));
    len = sqrt(length(param.phi));
    mk1 = (param.HState == 0) & Q & (param.RegionMask) & reshape(param.phi,len,len);

    % 2) Change from state 1 to 0
    mk0 = (param.HState == 1) & (param.H < param.H_thresh_off) & (param.RegionMask);

    % 3) Update
    % reset param values
%     collH(mk0) = 0;
%     collH(mk1) = param.H_thresh_on;
%     collA(mk0) = 0;
    param.H(mk0) = 0;
    param.H(mk1) = param.H_thresh_on;
    param.A(mk0) = 0;
    
    % reset states
    param.HState(mk0) = 0;
    param.HState(mk1) = 1;
    
%     param.A = reshape(NuNv\reshape(collA,sz(1)*sz(2),1),sz);
%     param.H = reshape(NuNv\reshape(collH,sz(1)*sz(2),1),sz);
    
    subplot(1,2,1);
    imagesc(param.A);
    colorbar;
    axis square;
    title('A');
    subplot(1,2,2);
    imagesc(param.H);
    colorbar;
    axis square;
    title('H');    
    drawnow;
    
    %% Get and write the frame
%     frame = getframe(1); % Grab frame
%     writeVideo(writer,frame); % Write frame to video

end

% close(writer); % Close videoWriter
% fprintf('Done!\n'); % Notify user when recording is complete

%% Functions
function dydt = evolvesystem(~,y,param)
% 0) unpack y = [A..H..]
L = length(y)/2;

sz_square = [sqrt(L) sqrt(L)];
A = reshape(y(1:L),sz_square);
H = reshape(y(L+1:2*L),sz_square);
phi_cp = reshape(param.phi,sz_square);

% 1) compute Short range Activator kernel
kSA = conv2((H.*phi_cp),param.G_SA,'same');

% 2) compute cytosolic hem pool
H_C = (param.HTotal - sum(sum(H)))/param.HTotal;

% 3) update differentials
dA = param.c_A(1).*(H.*phi_cp).*(1-A).*param.RegionMask;
dH = param.c_H(1).*((1-A).*H_C.*kSA - param.c_H(2)*A.*(H.*phi_cp)).*param.RegionMask;

% 4) non-nucleated hem sites don't evolve A or H
dA = dA.*param.HState;
dH = dH.*param.HState;

% 6) repack dY = [dA..dH..]
dydt = [reshape(dA,L,1); reshape(dH,L,1)].*[param.phi; param.phi];
% dydt = [reshape(dA,L,1); reshape(dH,L,1)];
end

% --------------------------------------------------------------------------
function param = GetParam
% load('phiplot_05202021');
load('phiplot_05202021_100by100');
len = length(phi_plot);
phi = full(phi_plot);
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
param.dt = 1;
param.T_current = 0;
param.W = 51;
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
