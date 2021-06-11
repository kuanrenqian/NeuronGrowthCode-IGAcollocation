function [phi,conc,conct,theta,tempr] = nucleus_tubulin_multiple(Nx,Ny,seed,epsilonb, center)

format long;

phi = zeros(Nx,Ny);
conc = 0.75.*ones(Nx,Ny);
conct = zeros(Nx,Ny);

for i=1:Nx
    for j=1:Ny
        for k=1:size(center,1)
            if ((i-center(k,1))*(i-center(k,1))+(j-center(k,2))*(j-center(k,2)) < seed)
                r = sqrt((i-center(k,1))*(i-center(k,1))+(j-center(k,2))*(j-center(k,2)));
                phi(i,j) = 1.0;
                conc(i,j) = 0.5;
                conct(i,j) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
            end
        end
    end
end

rng('shuffle');

theta = rand(Nx,Ny);

tempr = zeros(Nx,Ny);

end %endfunction