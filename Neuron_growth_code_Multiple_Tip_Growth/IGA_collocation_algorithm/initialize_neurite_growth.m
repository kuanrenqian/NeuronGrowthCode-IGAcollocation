function [phi, conc] = initialize_neurite_growth(seed_radius, lenu, lenv)

seed = (seed_radius)^2;

ind_i = [];
ind_j = [];

phi_val = [];

for i=1:lenu
    for j=1:lenv
        if ((i-lenu/2)*(i-lenu/2)+(j-lenv/2)*(j-lenv/2) < seed)
            %% single neuron
            ind_i(end+1) = i;
            ind_j(end+1) = j;
            phi_val(end+1) = 1;            

            %% multi neuron
%             ind_i(end+1) = i+15;
%             ind_i(end+1) = i-15;
%             ind_j(end+1) = j+15;
%             ind_j(end+1) = j-15;
%             
%             phi_val(end+1) = 1;            
%             phi_val(end+1) = 1;                    
            
        end
    end
end

%% Creating sparse matrix
phi = sparse(ind_i,ind_j,phi_val,lenu,lenv);
phi = full(reshape(phi,lenu*lenv,1)); % reshpae phi for calculation
conc = phi*0.5; % conc is 0.5 where phi is 1