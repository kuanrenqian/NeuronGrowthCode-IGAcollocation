function parameters = setparameters_neurite(Nx,Ny)

rho(1,1) = 5;
rho(2,1) = 10;
rho(3,1) = 20;

maxlevel = 2;

nelemx = Nx;
nelemy = Ny;

pU = 3;
pV = 3;

parameters = struct('pU',pU,'pV',pV,'maxlevel',maxlevel,'nelemx',nelemx,'nelemy',nelemy,'rho',rho);
end