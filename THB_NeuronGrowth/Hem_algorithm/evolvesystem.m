function dydt = evolvesystem(~,y,param)
% 0) unpack y = [A..H..]
L = length(y)/2;

sz_square = [sqrt(L) sqrt(L)];
A = reshape(y(1:L),sz_square);
H = reshape(y(L+1:2*L),sz_square);
phi_cp = reshape(param.phi,sz_square);

% 1) compute Short range Activator kernel
kSA = conv2(full(H.*phi_cp),param.G_SA,'same');

% 2) compute cytosolic hem pool
H_C = (param.HTotal - sum(sum(H)))/param.HTotal;

% 3) update differentials
dA = param.c_A(1).*(H.*phi_cp).*(1-A).*param.RegionMask;
dH = param.c_H(1).*((1-A).*H_C.*kSA - param.c_H(2)*A.*(H.*phi_cp)).*param.RegionMask;

% 4) non-nucleated hem sites don't evolve A or H
dA = dA.*param.HState;
dH = dH.*param.HState;

% 6) repack dY = [dA..dH..]
dydt = [reshape(dA,L,1); reshape(dH,L,1)].*[param.phi; param.phi];
% dydt = [reshape(dA,L,1); reshape(dH,L,1)];
end