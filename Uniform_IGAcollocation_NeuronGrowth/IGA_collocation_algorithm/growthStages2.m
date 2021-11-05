function [theta_ori] = growthStages2(phi_plot)
%#codegen
    tip = sum_filter(full(phi_plot),0);
    regionalMaxima = imregionalmax(full(tip));
    [Max_y,Max_x] = find(regionalMaxima);
    size_Max = length(Max_x);
    [lenu,lenv] = size(phi_plot);
    [theta_ori] = theta_rotate(lenu,lenv,Max_x,Max_y,size_Max);
end
