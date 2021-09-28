function [output] = to0(input)
[w,h] = size(input);
for i = 1:w
    for j = 1:h
        if (input(i,j)< 1e-6)
            output(i,j) = 0;
        end
end
end