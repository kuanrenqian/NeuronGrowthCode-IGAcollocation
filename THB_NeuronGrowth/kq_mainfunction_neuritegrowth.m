clc;
clear variables;
close all;

% path to Aishwarya THB functions
addpath('./setparameters');
addpath('./thbspline');
addpath('./iterationloop_funcs');

% diary 'log_Neuron_Growth'
rngSeed = rng('shuffle');
save('./postprocessing/rngSeed','rngSeed');

% suppress griddata warning (plotting dup points warning)
% warn_id = 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId';

% suppress rankDeficient warning (rows of zero in matrix)
warn_id = 'MATLAB:rankDeficientMatrix'; 
warning('off',warn_id)

%% Phase Field Simulation Variable Initialization
% time stepping variables
% dtime = 5e-3;
dtime = 5e-5;
end_iter = 35000;

% tolerance for NR method
tol = 1e-3;

% neuron growth variables
aniso = 6;
kappa= 4;
alph = 0.9; % changing name to alph cause alpha is a function
pix=4.0*atan(1.0);
gamma = 15.0;
tau = 0.3;
M_phi = 60;
M_theta = 0.5*M_phi;
s_coeff = 0.007;

delta = 0.1;
epsilonb = 0.04;

% Tubulin parameters
alpha_t = 0.001;
beta_t = 0.001;
Diff = 4;
% source_coeff = 0.012;
source_coeff = 0.05;

% Seed size
seed_radius = 5;

% Expanding domain parameters
BC_tol = 10;
expd_coef = 1.2;

disp('Simulation parameters initialization done!');

%% Domain Setup
Nx = 20;
Ny = 20;
dx = 1/Nx;
dy = 1/Ny;

parameters = setparameters_neurite(Nx,Ny);

% this initialization is purely for local refinement, actual phi is
% initialized in iteration_loop
[phi,conct] = kqInitializeNeuriteGrowth(seed_radius,Nx);

%%
tic
kq_MultipleResolution2D_neurite
toc