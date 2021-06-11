function [phi,conc,conct,theta,tempr] = nucleus_tubulin(Nx,Ny,seed,epsilonb)

format long;

%for i¼1:Nx
%for j¼1:Ny

%phi(i,j) ¼ 0.0;
%tempr(i,j) ¼ 0.0;

%end
%end

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

%rng('shuffle'); 

theta = rand(Nx,Ny);

tempr = zeros(Nx,Ny);

end %endfunction