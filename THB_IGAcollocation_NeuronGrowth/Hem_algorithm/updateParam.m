function param = updateParam(param,phi)
len = length(phi);
phi = full(phi);
p = phi;
for i=1:len
    if(phi(i)>0.5)
        p(i) = 1;
    else
        p(i) = 0;
    end
end
param.RegionMask = p;
param.phi = reshape(phi,len,1);
end
