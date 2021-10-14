function [epsilon_cp, epsilon_deriv_cp, aap_cp, P_dy, P_dx] = kqGetEpsilonAndAap(epsilonb,delta,phi_cp,theta,NuNv,NuN1v,N1uNv)
% kqGet_epsilon_and_aap(phi,theta,NuN1v,N1uNv) 
% This function calculates epsilon and aap (a*a') based on phi, theta,
% NuN1v, and N1uNv.
% kqGet_epsilon_and_aap() output epsilon aap in the format of [epsilon,aap]

aniso = 6;

P_dy = NuN1v*phi_cp;
P_dx = N1uNv*phi_cp;

atheta =(full(atan2(P_dy,P_dx)));

epsilon = epsilonb.*(1.0+delta*cos(aniso*(atheta-theta)));
epsilon_deriv = -epsilonb.*(aniso*delta*sin(aniso.*(atheta-theta)));
aap = epsilon.*epsilon_deriv;

% potential issue here, needs to be solved.
% epsilon_cp = NuNv\epsilon;
% epsilon_deriv_cp = NuNv\epsilon_deriv;
% aap_cp = NuNv\aap;

epsilon_cp = epsilon;
epsilon_deriv_cp = epsilon_deriv;
aap_cp = aap;
