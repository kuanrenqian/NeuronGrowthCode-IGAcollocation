clc
clear all
close all

path = pwd;
addpath(strcat(path,'\setparameters'));
addpath(strcat(path,'\thbspline'));
addpath(strcat(path,'\iterationloop'));
addpath(strcat(path,'\postprocessing'));

parameters = setparameters_neurite_NN2();

%% Read Geometry, Initialize
Nx = 350;
Ny = 350;
dx = 1;
dy = 1;

[laplacian_mat] = laplacian(Nx,Ny,dx,dy);
NxNy=(Nx)*(Ny);

%% Time integration parameters:
dtime = 5.0e-3;
dtime_theta = 1.e-3;
dtime_conc = 5.e-3;
dtime_conct = 5.e-3;

%% Material specific parameters:
mu = 1.0;
kappa= 1.8;
delta=0.5;
aniso= 2.0;
alpha = 0.9;
gamma = 50.0;
teq = 1.0;
theta0 = 0.2;
seed = (20*dx)^2;
gamma_w = 1;
delta_w = 6*dx;
lambda = 1.0e-16;
M = 10;
b = 88.9;
k_conc = 0.5;
D_cell = 10;
D_liquid = 1;
nstep = 200;
%% Calculate constant values
pix=4.0*atan(1.0);
abar = sqrt(3*delta_w*gamma_w/b);
epsilonb = (abar/(1+delta));
tau = epsilonb*dx;
W = (6*gamma_w*b)/delta_w;
M_phi = 50;
M_theta = 0.5*M_phi ;
epsilon0 = 0.012;

%% Tubulin parameters
alpha_t = 0.001;
beta_t = 0.001;
Diff = 1e2;
source_coeff = epsilon0;
iter_0 = 200;

%% Initialization of variables
center = [Nx/4,Ny/4;Nx/4,3*Ny/4;Nx/2,Ny/2;3*Nx/4,Ny/4;3*Nx/4,3*Ny/4];
flag_center = zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        if(i<=Ny/3 && j<=Nx/2)
            flag_center(i,j) = 0.1;
        elseif(i<=Ny/3 && j>=Nx/2)
            flag_center(i,j) = 0.2;
        elseif(i>=Ny/3 && i<=2*Ny/3)
            flag_center(i,j) = 0.3;
        elseif(i>=2*Ny/3 && j<=Nx/2)
            flag_center(i,j) = 0.4;
        elseif(i>=2*Ny/3 && j>=Nx/2)
            flag_center(i,j) = 0.5;
        end
    end
end
      
dist_center = sqrt((center(1,1)-center(2,1))^2+(center(1,2)-center(2,2))^2)

[phi,conc,conc_t,theta,tempr] = nucleus_tubulin_multiple(Nx,Ny,seed,epsilonb,center);
phi_original = phi;
phi3 = zeros(Nx,Ny);
flag_1 = zeros(Nx,Ny);

theta_original = theta;
save('theta_original_NN2.mat','theta_original');

% theta1 = load('theta_original_NN2.mat');
% theta = theta1.theta_original;

figure
imagesc(phi+flag_center)

MultipleResolution2D_neurite_NN2
