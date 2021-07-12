function [matdx,matdy] = gradient_mat_imf(mat,dx,dy)
%--
%-- this function uses imfilter to compute the discrete gradients
%--

% Prewitt filters
filter_x = [-1, 0, 1; -1, 0, 1; -1, 0, 1]/(6*dx);
filter_y = [-1, -1, -1; 0, 0, 0; 1, 1, 1]/(6*dy);

matdx = imfilter(mat, filter_x, 'circular');
matdy = imfilter(mat, filter_y, 'circular');


end
