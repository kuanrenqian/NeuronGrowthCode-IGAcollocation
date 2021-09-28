function [var_plot] = kqTHBplot(input)
[M,~] = size(input);
N = sqrt(M);
var_plot = reshape(input,N,N);
imagesc(var_plot);
colorbar;