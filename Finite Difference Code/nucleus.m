function [phi,tempr,theta, flag_theta] = nucleus_old(Nx,Ny,seed)

phi = zeros(Nx,Ny);
tempr = zeros(Nx,Ny);
theta1 = rand(Nx,Ny);
flag_theta = zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        if ((i-Nx/2)*(i-Nx/2)+(j-Ny/2)*(j-Ny/2) < seed)
            phi(i,j) = 1.0;
            flag_theta(i,j) = 1.0;
        end
    end
end
theta = theta1;
