clc
clear all
close all

path = pwd;
addpath(strcat(path,'\setparameters'));
addpath(strcat(path,'\thbspline'));
addpath(strcat(path,'\iterationloop'));
addpath(strcat(path,'\postprocessing'));

parameters = setparameters_neurite_2n();

%% Read Geometry, Initialize
Nx = 200;
Ny = 200;
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
kappa= 1.8;
delta=0.5;
aniso= 2.0;
alpha = 0.9;
gamma = 15.0;
seed = (40*dx)^2;
k_conc = 0.5;
D_cell = 10;
D_liquid = 1;
nstep = 600;

%% Calculate constant values
pix=4.0*atan(1.0);
abar = 0.45;
epsilonb = (abar/(1+delta));
tau = epsilonb;
M_phi = 50;
M_theta = 0.5*M_phi ;

%% Tubulin parameters
alpha_t = 0.001;
beta_t = 0.001;
Diff = 1e2;
source_coeff = 1.0;
iter_0 = 500;

%% Initialization of variables
center = [Nx/2,Ny/2];
[phi,conc,conc_t,~,tempr] = nucleus_tubulin_multiple(Nx,Ny,seed,epsilonb,center);
phi_original = phi;
phi3 = zeros(Nx,Ny);
flag_1 = zeros(Nx,Ny);
%theta_original = theta;
%save('theta_original_2n_2.mat','theta_original');

theta1 = load('theta_original_nc_3.mat');
theta = theta1.theta_original;

figure
imagesc(phi)

MultipleResolution2D_neurite_2n
