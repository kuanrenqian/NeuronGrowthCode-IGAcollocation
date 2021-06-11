function parameters = setparameters_neurite_NN2()

rho(1,1) = 5;
rho(2,1) = 2;
rho(3,1) = 5;

maxlevel = 3;

nelemx = 20;
nelemy = 20;

pU = 3;
pV = 3;

parameters = struct('pU',pU,'pV',pV,'maxlevel',maxlevel,'nelemx',nelemx,'nelemy',nelemy,'rho',rho);
end