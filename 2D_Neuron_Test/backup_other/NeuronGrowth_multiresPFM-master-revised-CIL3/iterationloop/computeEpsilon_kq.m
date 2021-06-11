function [epsilon,epsilon_deriv,dist1] = computeEpsilon(phi,phidx,phidy,theta,epsilonb,delta,aniso,Nx,Ny,epsilon,epsilon_deriv,dist1)
%#codegen

% Variable initialization moved to Multiresolution
% epsilon = zeros(Nx,Ny);
% epsilon_deriv = zeros(Nx,Ny);
% dist1 = zeros(Nx,Ny);

parfor i =1:Nx
    for j = 1:Ny

        %-- calculate angle:
        atheta =atan2(phidy(i,j),phidx(i,j));
        xtheta = 2*pi*theta(i,j);

        if(pi/4<=(atheta-xtheta)<=3*pi/4)
            epsilon(i,j) = epsilonb*(1.0+delta*cos(aniso*(atheta-xtheta)));
            epsilon_deriv(i,j) = -epsilonb*aniso*delta*sin(aniso.*(atheta-xtheta));

        elseif(-pi/4<(atheta-xtheta)<pi/4)
            epsilon(i,j) = (epsilon1/cos(pi/4)).*(cos((atheta-xtheta)));
            epsilon_deriv(i,j) = -(epsilon1/cos(pi/4)).*sin((atheta-xtheta));

        elseif(5*pi/4<=(atheta-xtheta)<=7*pi/4)
            epsilon(i,j) = epsilonb*(1.0+delta*cos(aniso*(atheta-xtheta)));
            epsilon_deriv(i,j) = -epsilonb*aniso*delta*sin(aniso.*(atheta-xtheta));

        elseif(3*pi/4<(atheta-xtheta)<5*pi/4)
            epsilon(i,j) = (epsilon1/cos(pi/4)).*(cos((atheta-xtheta)));
            epsilon_deriv(i,j) = -(epsilon1/cos(pi/4)).*sin((atheta-xtheta));
        end

        %evaluate the distance of neurites
        if(phi(i,j)>=0.5)
            dist1(i,j) = sqrt((i-Nx/2)^2+(j-Ny/2)^2);
        end
    end
end

end
