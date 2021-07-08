function [max_x,max_y,dist] = find_tip(phi_full,lenu,lenv)

dist= zeros(lenu,lenv);
for i = 1:lenu
    for j = 1:lenv
        if(phi_full(i,j)>0.85)
            dist(i,j) = sqrt((i-lenu/2)^2+(j-lenv/2)^2);
        end
    end
end

dist = reshape(dist,lenu*lenv,1);
[max_dist,max_index] = max(dist);
max_x = ceil(max_index/lenu);
max_y = rem(max_index,lenu);