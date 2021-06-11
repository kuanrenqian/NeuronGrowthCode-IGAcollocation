function parameters = setparameters_neurite()

rho(1,1) = 1;
rho(2,1) = 1;
rho(3,1) = 3;

maxlevel = 3;

nelemx = 80;
nelemy = 80;

pU = 3;
pV = 3;

parameters = struct('pU',pU,'pV',pV,'maxlevel',maxlevel,'nelemx',nelemx,'nelemy',nelemy,'rho',rho);
end