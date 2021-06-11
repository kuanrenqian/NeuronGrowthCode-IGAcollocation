clc
clear all
close all

path = pwd;
addpath(strcat(path,'\setparameters'));
addpath(strcat(path,'\thbspline'));
addpath(strcat(path,'\iterationloop'));
addpath(strcat(path,'\postprocessing'));
addpath(strcat(path,'\CIL_data'));
parameters = setparameters_neurite_CIL3();

%% Read Geometry, Initialize
Nx = 700;
Ny = 700;
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
seed = (28*dx)^2;
k_conc = 0.5;
D_cell = 10;
D_liquid = 1;
nstep = 400;

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
iter_0 = 200;

%% Initialization of variables
img = imread('8786_orig.tif');
img = img(:,:,1);
x1 = [465,226];
x2 = [725,772];
dist = sqrt((x1(1,1)-x2(1,1))^2 + (x1(1,2)-x2(1,2))^2);
distx = x1(1,1)-x2(1,1);
disty = x1(1,2)-x2(1,2);
figure
imagesc(img)
colormap gray

x1 = 65;
y1 = 250;
center = [x1, y1; x1+546, y1+260];
flag_center = zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        if(i<=Nx/2)
            flag_center(i,j) = 1;
        end
    end
end
      
dist_center = sqrt((center(1,1)-center(2,1))^2+(center(1,2)-center(2,2))^2)
dist_y = abs(center(1,1)-center(2,1))
dist_x = abs(center(1,2)-center(2,2))

save('center_CIL3.mat','center');

[phi,conc,conc_t,theta,tempr] = nucleus_tubulin_multiple(Nx,Ny,seed,epsilonb,center);
phi_original = phi;
phi3 = zeros(Nx,Ny);
flag_1 = zeros(Nx,Ny);
theta_original = theta;

save('theta_original_CIL3.mat','theta_original');

% theta1 = load('theta_original_CIL3.mat');
% theta = theta1.theta_original;
theta_original = theta;

delta_theta = pi/12;
for i =1:Nx
    for j=1:Ny
        theta2 =atan2(i-center(2,1),j-center(2,2));
        if(sqrt((i-center(1,1))^2+(j-center(1,2))^2)<=1.5*sqrt(seed))
            theta1 =atan2(i-center(1,1),j-center(1,2));
            if(theta1 >= pi/5-delta_theta && theta1 <= pi/5+delta_theta)
            elseif(theta1 >= 11*pi/18-delta_theta && theta1 <= 11*pi/18+delta_theta)  
            elseif(theta1 >= -pi/3-delta_theta && theta1 <= -pi/3+delta_theta)
            elseif(theta1 >= -13*pi/18-delta_theta && theta1 <= -13*pi/18+delta_theta)
            else
                theta(i,j) = 0;
            end
        end
        if(sqrt((i-center(2,1))^2+(j-center(2,2))^2)<=1.5*sqrt(seed))
            theta2 =atan2(i-center(2,1),j-center(2,2));
            if(theta2 >= pi/3-delta_theta && theta2 <= pi/3+delta_theta)
            elseif(theta2 >= 20*pi/36-delta_theta && theta2 <= 20*pi/36+delta_theta)  
            elseif(theta2 >= -4*pi/9-delta_theta && theta2 <= -4*pi/9+delta_theta) 
            elseif(theta2 >= -19*pi/24-delta_theta && theta2 <= -19*pi/24+delta_theta)     
            else
                theta(i,j) = 0;
            end
        end
    end
end


figure
imagesc(phi+flag_center)
figure
imagesc(theta)

MultipleResolution2D_neurite_CIL3
