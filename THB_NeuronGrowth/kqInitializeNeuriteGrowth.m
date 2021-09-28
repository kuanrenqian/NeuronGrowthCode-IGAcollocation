function [phi,conct] = kqInitializeNeuriteGrowth(seed_radius, coll_sz)

seed = (seed_radius)^2;

ind_i = [];
ind_j = [];

phi_val = [];
conct_val = [];

for i=1:coll_sz
    for j=1:coll_sz
        if ((i-coll_sz/2)*(i-coll_sz/2)+(j-coll_sz/2)*(j-coll_sz/2) < seed)
            r = sqrt((i-coll_sz/2)*(i-coll_sz/2)+(j-coll_sz/2)*(j-coll_sz/2));

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
            
            
            %% 3 neuron
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
%             ind_i(end+1) = i-45;
%             ind_i(end+1) = i-45;
%             ind_i(end+1) = i+45;
%             ind_i(end+1) = i+45;
% 
%             ind_j(end+1) = j-45;
%             ind_j(end+1) = j+45;
%             ind_j(end+1) = j+45;  
%             ind_j(end+1) = j-45;  
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
phi = sparse(ind_i,ind_j,phi_val,coll_sz,coll_sz);
phi = full(reshape(phi,coll_sz*coll_sz,1)); % reshpae phi for calculation
conct = sparse(ind_i,ind_j,conct_val,coll_sz,coll_sz);
conct = full(reshape(phi,coll_sz*coll_sz,1)); % reshpae phi for calculation
