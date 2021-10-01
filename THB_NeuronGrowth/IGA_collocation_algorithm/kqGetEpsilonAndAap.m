function [epsilon_cp, epsilon_deriv_cp, aap_cp, P_dy, P_dx] = kqGetEpsilonAndAap(epsilonb,delta,phi_cp,theta_cp,NuNv,NuN1v,N1uNv)
% kqGet_epsilon_and_aap(phi,theta,NuN1v,N1uNv) 
% This function calculates epsilon and aap (a*a') based on phi, theta,
% NuN1v, and N1uNv.
% kqGet_epsilon_and_aap() output epsilon aap in the format of [epsilon,aap]

aniso = 6;

theta = NuNv*theta_cp;
% phi_cp = NuNv\phi_cp;
P_dy = NuN1v*phi_cp;
P_dx = N1uNv*phi_cp;

% sz = sqrt(length(P_dx));
% 
% P_dy = reshape(P_dy,sz*sz,1);
% P_dx = reshape(P_dx,sz*sz,1);

atheta =(full(atan2(P_dy,P_dx)));
% atheta_cp = NuNv\atheta;

% xtheta = theta_cp;
% epsilon_cp = epsilonb.*(1.0+delta*cos(aniso*(atheta-xtheta)));
% epsilon_deriv_cp = -epsilonb.*(aniso*delta*sin(aniso.*(atheta_cp-xtheta)));
% epsilon_cp = epsilonb.*(1.0+delta*cos(aniso*(atheta_cp-theta_cp)));
% epsilon_deriv_cp = -epsilonb.*(aniso*delta*sin(aniso.*(atheta_cp-theta_cp)));
% aap_cp = epsilon_cp.*epsilon_deriv_cp;

epsilon = epsilonb.*(1.0+delta*cos(aniso*(atheta-theta)));
epsilon_deriv = -epsilonb.*(aniso*delta*sin(aniso.*(atheta-theta)));
aap = epsilon.*epsilon_deriv;

epsilon_cp = NuNv\epsilon;
epsilon_deriv_cp = NuNv\epsilon_deriv;
aap_cp = NuNv\aap;
