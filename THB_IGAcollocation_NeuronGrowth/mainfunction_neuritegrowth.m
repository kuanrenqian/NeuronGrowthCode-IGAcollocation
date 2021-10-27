clc;
clear all;
close all;

addpath('./setparameters');
addpath('./thbspline');
addpath('./iterationloop_funcs');

% suppress griddata warning (plotting dup points warning)
warn_id = 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId';
warning('off',warn_id)

Nx = 20;
Ny = 20;
dx = 1/Nx;
dy = 1/Ny;

rngSeed = rng('shuffle');

parameters = setparameters_neurite(Nx,Ny);

[phi,conct] = kqInitializeNeuriteGrowth(5,Nx);
phi = reshape(phi,Nx,Ny);

tic
MultipleResolution2D_neurite
toc