function [phi, conc, conct] = initialize_neurite_growth(seed_radius, lenu, lenv, dx)
phi = zeros([lenu,lenv]); % with ghost nodes (technically Nx by Ny, but becomes size_collpts after adding ghost points)
conc = 0.75.*ones([lenu,lenv]);
conct = zeros(lenu,lenv);
seed = (seed_radius*dx)^2;
for i=1:lenu
    for j=1:lenv
        if ((i-lenu/2)*(i-lenu/2)+(j-lenv/2)*(j-lenv/2) < seed)
            r = sqrt((i-lenu/2)*(i-lenu/2)+(j-lenv/2)*(j-lenv/2));
%             phi(i+15,j+15) = 1.0;
%             phi(i-15,j-15) = 1.0;
%             conc(i+15,j+15) = 0.5;
%             conc(i-15,j-15) = 0.5;
            
            phi(i,j) = 1.0;
            conc(i,j) = 0.5;
            conct(i,j) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
        end
    end
end