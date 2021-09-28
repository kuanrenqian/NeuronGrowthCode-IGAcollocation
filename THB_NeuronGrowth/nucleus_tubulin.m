function [phi,conc,conct,theta,tempr] = nucleus_tubulin(Nx,Ny,seed,epsilonb)

format long;

phi = zeros(Nx,Ny);
conc = 0.75.*ones(Nx,Ny);
conct = zeros(Nx,Ny);

for i=1:Nx
    for j=1:Ny
        if ((i-Nx/2)*(i-Nx/2)+(j-Ny/2)*(j-Ny/2) < seed)
            r = sqrt((i-Nx/2)*(i-Nx/2)+(j-Ny/2)*(j-Ny/2));
            phi(i,j) = 1.0;
            conc(i,j) = 0.5;
            conct(i,j) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
        end
    end
end

rng('shuffle'); 

theta = randn(Nx,Ny);

tempr = zeros(Nx,Ny);

end %endfunction