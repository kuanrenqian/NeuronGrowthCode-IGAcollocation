function [grad] =laplacian_imf(m,dx,dy)
% Compute the discrete laplacian of the matrix m using imfilter

% 5-point stencil
lap_filter = [0, 1, 0; 1, -4, 1; 0, 1, 0]/(dx*dy);

% 9-point stencil
%lap_filter = [0.25, 0.5, 0.25; 0.5, -3, 0.5; 0.25, 0.5, 0.25]/(dx*dy);

grad = imfilter(m, lap_filter, 'replicate');

end 
