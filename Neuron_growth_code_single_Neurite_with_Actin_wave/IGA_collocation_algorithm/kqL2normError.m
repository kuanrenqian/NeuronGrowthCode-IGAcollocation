function L2error = kqL2normError(Analytical,Calculated)
% L2norm_error calculates L2 norm error for the input matrices

size_a = size(Analytical);
size_c = size(Calculated);

if(size_a~=size_c)
    error('Input matrix dimention mismatch!');
end

L2error = (Analytical-Calculated).^2;
L2error = sum(sum(L2error));
L2error = sqrt(L2error);
L2error = L2error/(size_a(1)*(size_a(2)));