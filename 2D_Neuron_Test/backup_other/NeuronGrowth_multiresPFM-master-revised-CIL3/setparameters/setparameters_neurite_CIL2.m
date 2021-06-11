function parameters = setparameters_neurite_CIL2()

rho(1,1) = 1;
rho(2,1) = 1;
rho(3,1) = 2;

maxlevel = 3;

nelemx = 120;
nelemy = 120;

pU = 3;
pV = 3;

parameters = struct('pU',pU,'pV',pV,'maxlevel',maxlevel,'nelemx',nelemx,'nelemy',nelemy,'rho',rho);

end