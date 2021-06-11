clc
clear all
close all
format long

addpath('utils');

% dx2 = 0.075;
% dy2 = 0.075;
dx2 = 0.03;
dy2 = 0.03;
nx2 = 100;
ny2 = 100;
dV2 = dx2*dy2;

Cle = 0.394;
Cse = 0.5413;
k = Cse/Cle;
initConc = Cle + 0.50*(Cse-Cle);
% delta=0.5;
aniso= 6.0;
% abar = 0.45;
% epsilonb = 0.5;
% epsilonb = 0.01;
% epsilonb = 1;

delta = 0.02;
epsilonb = 0.01;
tau = 0.0003;

x2 = 0:dx2:nx2*dx2;
y2 = 0:dy2:ny2*dy2;
[X,Y] = meshgrid(x2,y2);
phi = zeros(size(X,1),size(X,2));
conc = initConc*ones(size(X,1),size(X,2));
for i = 1:nx2
    for j = 1:ny2
%         if(((X(j,i) - (nx2/2-5)*dx2)^2 + (Y(j,i) - (ny2/2-5)*dy2)^2)<(15*dx2)^2)
%             phi(j,i) = 1;
%         end
%         if(((X(j,i) - (nx2/2+12)*dx2)^2 + (Y(j,i) - (ny2/2+12)*dy2)^2)<(8*dx2)^2)
%             phi(j,i) = 1;
%         end
        if(((X(j,i) - (nx2/2)*dx2)^2 + (Y(j,i) - (ny2/2)*dy2)^2)<(10*dx2)^2)
            phi(j,i) = 1;
        end
    end
end

tempr = zeros(size(X,1),size(X,2));

figure(1);
imagesc(phi)
colormap jet
colorbar;
title('initial phi');

%% Look Up Table for interpolation of Cs and Cl
LUTnp = 1002;
LUTnc = 1002;
LUTdp = 1./(LUTnp-2);
LUTdc = 1./(LUTnp-2);
flag1 = 1;

if(flag1==0)
    pureCs = zeros(LUTnp, LUTnc);
    pureCl = zeros(LUTnp, LUTnc);
    
    for i1 = 1:LUTnp
        for  j1 = 1:LUTnc
            LUTp = (i1-1)*LUTdp;
            LUTc = (j1-1)*LUTdc;
            [pureCs(i1,j1),pureCl(i1,j1)] = commonTangentConc(LUTp, LUTc, LUTc, LUTc);
        end
    end
else
    load 'pureCs.mat';
    load 'pureCl.mat';
end

figure(2);
subplot(2,1,1);
imagesc(pureCs);

subplot(2,1,2);
imagesc(pureCl);

LUTp  = linspace(-LUTdp, 1+LUTdp, LUTnp);
LUTc  = linspace(-LUTdc, 1+LUTdc, LUTnc);
[LUX, LUY] = meshgrid(LUTc,LUTp);
%[cs, cl] = interpolateCsCl(pureCs,pureCl,LUX,LUY, phi, conc);

smoothpureCs = interp2(pureCs,LUX,LUY, 'spline');
smoothpureCl = interp2(pureCl,LUX,LUY, 'spline');

% figure(3);
% subplot(2,1,1);
% imagesc(smoothpureCs);
% colorbar;
% colormap jet
% title('smooth pure Cs');
% 
% subplot(2,1,2);
% imagesc(smoothpureCl);
% colorbar;
% colormap jet
% title('smooth pure Cl');

Cs = zeros(size(phi,1),size(phi,2));
Cl = zeros(size(phi,1),size(phi,2));
for i =1:size(phi,1)
    for j = 1:size(phi,2)
        [Cs(i,j),Cl(i,j)] = commonTangentConc(phi(i,j), conc(i,j), conc(i,j), conc(i,j));
    end
end
% 
% figure(4);
% subplot(2,1,1);
% imagesc(Cs)
% colorbar;
% colormap jet
% title('initial Cs');
% 
% subplot(2,1,2);
% imagesc(Cl)
% colormap jet
% colorbar
% title('initial Cl');

%% Now start the solver for phase field
dtime = 1e-4;
nstep = 5000;
epsSq = 1.25;
halfWidth =5.0;

for istep = 1:nstep

    istep

    phiold = phi;
    concold = conc;
    Csold = Cs;
    Clold = Cl;
    
    temprold = tempr;

    alpha = 0.9;
    pi = 3.1415926;
    teq = 1;
    gamma = 10;
    % E - Energy
    m = alpha/pi * atan(gamma*(teq-temprold)); 

    %-- calculate angle:
    [phidy, phidx] = gradient_mat_imf(phiold,dx2, dy2);
    atheta = atan2(phidy,phidx);
    xtheta = 0.02;
        
    epsilon = epsilonb*(1.0+delta*cos(aniso*(atheta-xtheta)));
    epsilon_deriv = -epsilonb*aniso*delta*sin(aniso.*(atheta-xtheta));
    
    %calculate first term in phi equation:
    dummyx = epsilon.*epsilon_deriv.*phidx;
    [term1,~] = gradient_mat_imf(dummyx,dx2,dy2);
    
    %calculate second term in phi equation:
    dummyy =-epsilon.*epsilon_deriv.*phidy;
    [~,term2] =gradient_mat_imf(dummyy,dx2,dy2);
    
    %Cll = zeros(size(Cl,1),1);
%     omegaVal = omega(dx2);
    omegaVal = omega(dx2,epsilon);
    fLCL = fL(Cl);
    fSCS = fS(Cs);
    dfLdC = dfLdc(Cl);

    Gprime = linearGPrime(phiold);
    Gdoubleprime = linearGdoublePrime(phiold);
    Hprimee = linearHPrime(phiold).*tempr;
	Hdoubleprimee = linearHdoublePrime(phiold).*tempr;

%     Gp = 2.*phiold-6.*phiold.^2+4.*phiold.^3;
%     Gpp = 2 - 12.*phiold + 12.*phiold.^2;
%     Hp = 30.*phiold.^2 - 60.*phiold.^3 + 30.*phiold.^4;
%     Hpp = 60.*phiold - 180.*phiold.^2 + 120.*phiold.^3;
% 
%     Hp = -phiold.^3 + (1.5-m).*phiold.^2 + (m-0.5).*phiold;
%     Hpp = -3.*phiold.^2 + 2.*(1.5-m).*phiold + (m-0.5);
   
    diffusionPotential = fLCL - fSCS - dfLdC.*(Cl-Cs);

    mPhi = -omegaVal.*Gprime + Hprimee.*diffusionPotential;

    dmPhidPhi = -omegaVal.*Gdoubleprime + Hdoubleprimee.*diffusionPotential;

    S1 = dmPhidPhi.*phiold.*(1-phiold) + mPhi.*(1-2.*phiold);

    S0 = mPhi.*phiold.*(1-phiold) - S1.*phiold;

    lap_phi = laplacian_imf(phiold, dx2, dy2);
    
%     delta_phi = epsSq.*lap_phi + S0 + S1.*phiold;
%    delta_phi = epsilon.^2.*lap_phi + S0 + S1.*phiold;
    delta_phi = term1 +term2 + epsilon.^2.*lap_phi + S0 + S1.*phiold;
%     f_phi = -temprold.*Hp.*diffusionPotential+omegaVal.*Gp;
%     delta_phi = term1 +term2 + epsilon.^2.*lap_phi + phiold.*(1.0-phiold).*(phiold - 0.5 + m);
%     delta_phi = term1 +term2 + epsilon.^2.*lap_phi - f_phi;
%     delta_phi = term1 +term2 + epsilon.^2.*lap_phi + Hp.*diffusionPotential;

%     xi = -1 + 2.*rand(nx2+1,ny2+1);
%     pertub = 0.02;    
%     gphi = phiold.^2.*(1-phiold.^2);
%     
%     pertFun = 16.*xi.*gphi.*pertub;
% %     delta_phi = delta_phi + pertFun;

    phi = phiold + dtime.*delta_phi/tau;
%     phi = phiold + dtime.*delta_phi;
    
   %% concentration
   k = Cs./Cl;
    Qphi = evaluateQ(phiold, k);
    Dc = Qphi;
    Dp = Qphi.*Hprime(phiold).*(Cl-Cs);

    lap_conc = laplacian_imf(concold, dx2, dy2);
    diff_term = Dc.*lap_conc;
    [phidy, phidx] = gradient_mat_imf(phiold,dx2, dy2);

    [~, explicitsource_term_x] = gradient_mat_imf(Dp.*(phidx),dx2,dy2);
    [explicitsource_term_y, ~] = gradient_mat_imf(Dp.*(phidy),dx2,dy2);
    explicitsource_term = explicitsource_term_x + explicitsource_term_y;

    delta_conc = diff_term + explicitsource_term;

    conc = concold + dtime.*delta_conc;

    for i =1:size(phi,1)
        for j = 1:size(phi,2)
        [Cs(i,j),Cl(i,j)] = commonTangentConc(phi(i,j), conc(i,j), conc(i,j), conc(i,j));
        end
    end
    
    %% temperature evolve
    kappa = 2;
%     kk = 192.6;
%     rho = 2800;
%     Cp = 1130;
%     alph = kk/(rho*Cp);
%     L = 389*1000;
    lap_tempr = laplacian_imf(temprold, dx2, dy2);
    tempr = tempr + dtime*lap_tempr + kappa*(phi-phiold);
%     tempr = tempr + alph.*dtime*lap_tempr + L./Cp.*Hprime(phiold).*(phi-phiold);

    figure(7);
    subplot(4,2,1)
    imagesc(phi)
    colormap jet
    colorbar;
    title('current phi_v1');

    subplot(4,2,2)
    imagesc(tempr)
    colorbar;
    title('current tempr');

    subplot(4,2,3)
    imagesc( S0 + S1.*phiold)
    colorbar;
%     caxis([0,1]);
    title('current  S0 + S1.*phiold');

    subplot(4,2,4)
    imagesc(Cl)
    colorbar;
    caxis([0,1]);
    title('current Cl');
    
    subplot(4,2,5)
    imagesc(epsilon.^2)
    colorbar;
    title('epsilon.^2');
    
    subplot(4,2,6)
    imagesc(term1)
    colorbar;
    title('term1');
    
    subplot(4,2,7)
    imagesc(term2)
    colorbar;
    title('term2');
    
    subplot(4,2,8)
    imagesc(m)
    colorbar;
    title('E');
    
    drawnow
    
end