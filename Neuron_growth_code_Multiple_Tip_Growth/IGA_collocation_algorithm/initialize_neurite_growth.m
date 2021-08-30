function [phi,conct] = initialize_neurite_growth(seed_radius, lenu, lenv)

seed = (seed_radius)^2;

ind_i = [];
ind_j = [];

phi_val = [];
conct_val = [];

for i=1:lenu
    for j=1:lenv
        if ((i-lenu/2)*(i-lenu/2)+(j-lenv/2)*(j-lenv/2) < seed)
            r = sqrt((i-lenu/2)*(i-lenu/2)+(j-lenv/2)*(j-lenv/2));

            %% single neuron
            ind_i(end+1) = i;
            ind_j(end+1) = j;
            phi_val(end+1) = 1;            
            conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));

            %% 2 neurons
%             ind_i(end+1) = i+25;
%             ind_i(end+1) = i-25;

%             ind_j(end+1) = j+25;
%             ind_j(end+1) = j-25;
%             
%             phi_val(end+1) = 1;            
%             phi_val(end+1) = 1;                    
%             conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
%             conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));

            %% 3 neurons
%             ind_i(end+1) = i-30;
%             ind_i(end+1) = i+30;
%             ind_i(end+1) = i+30;
% 
%             ind_j(end+1) = j;
%             ind_j(end+1) = j+35;
%             ind_j(end+1) = j-35;  
%             
%             phi_val(end+1) = 1;            
%             phi_val(end+1) = 1;                    
%             phi_val(end+1) = 1;      
%             
%             conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
%             conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
%             conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));

            %% 4 neurons
%             ind_i(end+1) = i-35;
%             ind_i(end+1) = i-35;
%             ind_i(end+1) = i+35;
%             ind_i(end+1) = i+35;
% 
%             ind_j(end+1) = j-35;
%             ind_j(end+1) = j+35;
%             ind_j(end+1) = j+35;  
%             ind_j(end+1) = j-35;  
%             
%             phi_val(end+1) = 1;            
%             phi_val(end+1) = 1;                    
%             phi_val(end+1) = 1;      
%             phi_val(end+1) = 1;      
%             
%             conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
%             conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
%             conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
%             conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));

        end
    end
end

%% Creating sparse matrix
phi = sparse(ind_i,ind_j,phi_val,lenu,lenv);
phi = full(reshape(phi,lenu*lenv,1)); % reshpae phi for calculation
conct = sparse(ind_i,ind_j,conct_val,lenu,lenv);
conct = full(reshape(phi,lenu*lenv,1)); % reshpae phi for calculation
