clc;
clear all;
close all;

addpath('./setparameters');
addpath('./thbspline');
addpath('./iterationloop_funcs');

Nx = 20;
Ny = 20;
dx = 1/Nx;
dy = 1/Ny;

rngSeed = rng('shuffle');

parameters = setparameters_neurite(Nx,Ny);

[phi,conct] = kqInitializeNeuriteGrowth(5,22);
phi = reshape(phi,22,22);

tic
MultipleResolution2D_neurite
toc