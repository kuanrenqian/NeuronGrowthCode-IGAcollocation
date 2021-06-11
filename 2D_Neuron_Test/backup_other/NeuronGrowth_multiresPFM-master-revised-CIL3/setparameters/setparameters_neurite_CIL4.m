function parameters = setparameters_neurite_CIL4()

rho(1,1) = 1;
rho(2,1) = 4;
rho(3,1) = 6;

maxlevel = 3;

nelemx = 200;
nelemy = 200;

pU = 3;
pV = 3;

parameters = struct('pU',pU,'pV',pV,'maxlevel',maxlevel,'nelemx',nelemx,'nelemy',nelemy,'rho',rho);
end