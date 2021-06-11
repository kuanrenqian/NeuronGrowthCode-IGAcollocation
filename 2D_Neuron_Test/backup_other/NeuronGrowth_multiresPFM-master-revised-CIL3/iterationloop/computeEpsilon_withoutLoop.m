function [epsilon,epsilon_deriv,dist1] = computeEpsilon_withoutLoop(phi,phidx,phidy,theta,epsilonb,delta,aniso,Nx,Ny,epsilon1,epsilon,epsilon_deriv,dist1)

atheta_1 =atan2(phidy,phidx);
xtheta_1= 2*pi*theta;

index1 = (pi/4<=(atheta_1-xtheta_1)<=3*pi/4);
epsilon(index1) = epsilonb*(1.0+delta*cos(aniso*(atheta_1(index1)-xtheta_1(index1))));
epsilon_deriv(index1) = -epsilonb*aniso*delta*sin(aniso.*(atheta_1(index1)-xtheta_1(index1)));

index2 = (-pi/4<(atheta_1-xtheta_1)<pi/4)&(index1~=1);
epsilon(index2) = (epsilon1/cos(pi/4)).*(cos((atheta_1(index2)-xtheta_1(index2))));
epsilon_deriv(index2) = -(epsilon1/cos(pi/4)).*sin((atheta_1(index2)-xtheta_1(index2)));

index3 = (5*pi/4<=(atheta_1-xtheta_1)<=7*pi/4)&(index2~=1)&(index1~=1);
epsilon(index3) = epsilonb*(1.0+delta*cos(aniso*(atheta_1(index3)-xtheta_1(index3))));
epsilon_deriv(index3) = -epsilonb*aniso*delta*sin(aniso.*(atheta_1(index3)-xtheta_1(index3)));

index4 = (3*pi/4<(atheta_1-xtheta_1)<5*pi/4)&(index3~=1)&(index2~=1)&(index1~=1);
epsilon(index4) = (epsilon1/cos(pi/4)).*(cos((atheta_1(index4)-xtheta_1(index4))));
epsilon_deriv(index4) = -(epsilon1/cos(pi/4)).*sin((atheta_1(index4)-xtheta_1(index4)));

%evaluate the distance of neurites
dist1(phi>=0.5) = sqrt((i-Nx/2)^2+(j-Ny/2)^2);

end