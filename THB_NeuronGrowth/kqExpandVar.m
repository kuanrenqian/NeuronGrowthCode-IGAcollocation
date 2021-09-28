function [output] = kqExpandVar(input,Nx,Ny,isRand)
    input = reshape(input,Nx,Ny);
    if isRand == 1
        temp = rand(2*Nx,2*Ny);
    else
        temp = zeros(2*Nx,2*Ny);
    end
    temp(floor((Nx*2 -Nx)/2)+1:Nx+floor((Nx*2-Nx)/2),floor((Ny* 2  -Ny)/2)+1:Ny+floor((Ny* 2  -Ny)/2)) = input;

    Nx = Nx*2;
    Ny = Ny*2;

    output = reshape(temp,Nx*Ny,1);
end

%% Backup code
%         phi = reshape(phi,Nx,Ny);
%         phi_new = zeros(2*Nx,2*Ny);
%         phi_new(floor((Nx*2 -Nx)/2)+1:Nx+floor((Nx*2-Nx)/2),floor((Ny* 2  -Ny)/2)+1:Ny+floor((Ny* 2  -Ny)/2)) = phi;
% 
%         Nx = Nx*2;
%         Ny = Ny*2;
%         NxNy = (Nx)*(Ny);
%         
%         phi = reshape(phi_new,Nx*Ny,1);