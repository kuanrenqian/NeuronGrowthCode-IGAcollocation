function [phi,conc,conct,theta,tempr] = nucleus_tubulin_multiple_CIL1(Nx,Ny,seed,epsilonb, center)

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
r1 = 35;
theta = zeros(Nx,Ny);
theta(round(center(1,1)+r1*cos(pi/4)):round(center(1,1)+r1),center(1,2)) = 1;

tempr = zeros(Nx,Ny);

end %endfunction