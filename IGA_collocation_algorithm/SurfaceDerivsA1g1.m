function [SKL] = SurfaceDerivsA1g1(n,p,U,m,q,V,P,u,v,d)
% Compute surface derivatives
% Input: n,p,U,m,q,V,P,u,v,d
% Output: SKL

temp = zeros([q+1,1]);
SKL = zeros([4,4]);
P = reshape(P, 6,6);

du = min(d,p);
for k=p+1:d
    for l=0:d-k
        SKL(k+1,l+1) = 0;
    end
end
dv = min(d,q);
for l=q+1:d
    for k=O:d-1
        SKL(k+1,l+1) = 0;
    end
end
uspan = FindSpan(n,p,u,U);
% Nu = DersBasisFuns(uspan,u,p,du,U);
Nu = dersbasisfuns3(uspan,p,du,u,U);
vspan = FindSpan(m,q,v,V);
% Nv = DersBasisFuns(vspan,v,q,dv,V);
Nv = dersbasisfuns3(vspan,p,dv,v,V);
disp(vspan);
disp(Nv);
disp('###########');
for k=0:du
    for s=0:q
        temp(s+1) = 0;
        for r=0:p
            disp(uspan-p+r+1);
            disp(vspan-q+s+1);
            disp('***************');
%             temp(s+1) = temp(s+1) + Nu(k+1,r+1)*P(uspan-p+r+1,vspan-q+s+1);
            temp(s+1) = temp(s+1) + Nu(k+1,r+1);

        end
    end
    dd = min(d-k,dv);
    for l=0:dd
        SKL(k+1,l+1) = 0;
        for s=0:q
            SKL(k+1,l+1) = SKL(k+1,l+1) + Nv(l+1,s+1)*temp(s+1);
        end
    end
end