% Actin Wave model on 2D plate
% 05/18/2021


%% CleanUp
close all;
clear;
clc;

%% Including Path
addpath('../IGA_collocation_algorithm');

disp('********************************************************************');
disp('2D Reaction Diffusion Gray Scott solver using IGA-Collocation');

%% Variable Initialization
disp('Initializing vairbles...');
% time stepping variables
t = 0;
dt = 0.0005;
% dt = 1;

s1 = 0.5;
s2 = 1;
k0 = 0.2;
n = 3;

gamma = 1;
A0 = 0.4;
delta = 1;
F0 = 0.4;
kn = 1;
ks = 0.25;
Tnpf =1;
e = 0.1;
L = 1;
Da = (10e-3)/3;
Di = (1e-1)/3;

% B-spline curve order (U,V direction)
p = 3;
q = 3;

% Control mesh size (knot vector)
Nx = 60;
Ny = 60;
dx = 1;
dy =1;
knotvectorU = [zeros([1,p]),0:Nx,ones([1,p])*Nx].';
knotvectorV = [zeros([1,p]),0:Ny,ones([1,p])*Ny].';
lenu = length(knotvectorU)-2*(p-1);
lenv = length(knotvectorV)-2*(p-1);

%% Constructing coef matrix (concise function version)
order_deriv = 2;
sprs = 1;
[NuNv,N1uNv,NuN1v,~,N2uNv,NuN2v,~,~,size_collpts,~] = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv,sprs);
lap = N2uNv + NuN2v;

%% Iterating variable initialization
disp('Variabel initialization...');
a = zeros(size_collpts);
i = zeros(size_collpts);
f = zeros(size_collpts);

for k = floor(size_collpts/2)-5:floor(size_collpts/2)+5
    for j = floor(size_collpts/2)-5:floor(size_collpts/2)+5
        a(k,j) = 0.8;
%         i(k,j) = 1;
    end
end

a = reshape(a,lenu*lenv,1);
i = reshape(i,lenu*lenv,1);
f = reshape(f,lenu*lenv,1);

a_initial = a;
i_initial = i;
f_initial = f;

a = NuNv\a;
i = NuNv\i;

%% Parameters for Actin wave in literature


%% ID for boundary location (suppress 4 edges)
% % id = 1 means there is bc
% bcid = zeros([lenu,lenv]);
% for i = 1:lenu
%     bcid(1,i) = 1;
%     bcid(lenu,i) = 1;
%     bcid(i,1) = 1;
%     bcid(i,lenv) = 1;
% end
% bcid = sparse(reshape(bcid,lenu*lenv,1));
% N = zeros(lenu);


%% Start video writer
% writer = VideoWriter('active_wave_05202021'); % Define VideoWriter object
% open(writer);

%% Set figure(1) size and position
figure(1);
set(gcf,'position',[700,100,800,400]);

%%
for iter = 1:5000
    fprintf('Progress: %.2d\n',iter);

    %% Grey-Scott
    A = NuNv*a;
    dAdx = N1uNv*a;
    dAdy = NuN1v*a;
    lapA = (lap*a);
    I = NuNv*i;
    dIdx = N1uNv*i;
    dIdy = NuN1v*i;
    lapI = (lap*i);    
    F = NuNv*f;
    dFdx = N1uNv*f;
    dFdy = NuN1v*f;
    lapF = (lap*f);    
    
    %% Coupled equation
    
    f_new = NuNv\(e.*(kn.*A-ks.*F)*dt+F);

%     F = NuNv*f_new;
    
    term_fI = k0+(gamma.*A.^n)./(A0.^n+A.^n);
    term_fI = term_fI.*I;
    term_fA = delta.*(s1+s2.*(F./(F0+F)));
    term_fA = term_fA.*A;
    ff = term_fI - term_fA;
    
    dA = ff + Da .* lapA;
    dI = -ff + Di .* lapI;
    
    a_new = NuNv\(A + dt*dA);
    i_new = NuNv\(I + dt*dI);

%     uRHS  = (NuNv*u + dt*duT./P_temp);
%     uLHS = NuNv;
%     [uLHS, uRHS] = StiffMatSetupBCID(uLHS, uRHS,bcid,u_initial);
% 
% 
%     vRHS  = (NuNv*v + dt*dvT./P_temp);
%     vLHS = NuNv;
%     [vLHS, vRHS] = StiffMatSetupBCID(vLHS, vRHS,bcid,v_initial);
% 
%     u_new = uLHS\uRHS;
%     v_new = vLHS\vRHS;
    
%% Plotting figures
    if(mod(iter,1) == 0 )
        a_plot = reshape(NuNv*a_new,size_collpts,size_collpts);
        subplot(1,2,1);
        imagesc(a_plot);
        title(sprintf(' a iteration = %.2d',iter));
        axis square;
        colorbar;
        
        subplot(1,2,2);
        imagesc(reshape(NuNv*i_new,lenu,lenv));
        title('i');
        axis square;
        colorbar;

    %% plot current iteration
        drawnow;
        
    %% Get and write the frame
%         frame = getframe(1); % Grab frame
%         writeVideo(writer,frame); % Write frame to video
        
    end

    %% iteration update
    % update variables in this iteration    
    a = a_new;
    i = i_new;
    f = f_new;


    t = t+dt;
    
end

% close(writer); % Close videoWriter
% fprintf('Done!\n'); % Notify user when recording is complete

disp('All simulations complete!\n');
disp('********************************************************************');
