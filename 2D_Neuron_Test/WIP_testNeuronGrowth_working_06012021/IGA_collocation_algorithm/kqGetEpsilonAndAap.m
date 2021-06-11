function [epsilon, epsilon_deriv, aap] = kqGetEpsilonAndAap(epsilonb,delta,phi,theta,NuNv,NuN1v,N1uNv)
% kqGet_epsilon_and_aap(phi,theta,NuN1v,N1uNv) 
% This function calculates epsilon and aap (a*a') based on phi, theta,
% NuN1v, and N1uNv.
% kqGet_epsilon_and_aap() output epsilon aap in the format of [epsilon,aap]

aniso = 6;

P_dy = NuN1v*phi;
P_dx = N1uNv*phi;

% atheta =smooth(full(atan2(phidy,phydx)));
atheta =(full(atan2(P_dy,P_dx)));
xtheta = 2*pi*(NuNv*theta);
% xtheta = 0.2;

epsilon = epsilonb*(1.0+delta*cos(aniso*(atheta-xtheta)));
epsilon_deriv = -epsilonb*aniso*delta*sin(aniso.*(atheta-xtheta));
aap = epsilon.*epsilon_deriv;

epsilon = NuNv\epsilon;
aap = NuNv\aap;

epsilon = sparse(epsilon);
epsilon_deriv = sparse(epsilon_deriv);
aap = sparse(aap);
