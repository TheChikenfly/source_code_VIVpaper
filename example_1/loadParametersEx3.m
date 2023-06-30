% This file contains the parameters fo the hydro_frame_force example
% For Validation Examples 1 2 3
% Time parameters
%
global massratio;
finalTime = 1; dt = 0.00005;
%--------
%
ILVIVBool=true;
staticBool = false; % To start from the static deformed configuration
% Plot parameters
%
spanPlot = finalTime/dt/60; 
lw = 2.0 ; ms = 11 ; plotfontsize = 20 ;
%--------
% Fluid properties
%
rhoFluid = 1000; nuFluid = 1; 
% absolute mean velocity of water at x=[0 0 0] and time 0 
nameFuncVel = 'windVelUniform3'; 
va_vector = feval(nameFuncVel, [0, 0, 0], 0);
va        = norm(va_vector);
% extract lift and drag which are constant (betaRel and Re doesen´t matter)
nameLiftFunc = 'liftCoefWOM'; 
nameDragFunc = 'dragCoefWOM'; 
cL0 = feval(nameLiftFunc, 0, 1); cD0 = feval(nameDragFunc, 0, 1); 
cDi = 0.1;
%--------
%
% Material and geometric properties
%
% material
rho = rhoFluid; E = 5e8 ; nu = .3; 
massratio = rho/rhoFluid ;
% geometry
%
l = 0.2 ; d = 0.01; I = pi * d^4 / 64 ;
% stiffness
ku = 3*E*I/(l^3); 
% nodal lumped mass
m  = rho*pi*(d^2)*l/(4*2) ;
%Try to respect paper
%m  = (rho*pi*(d^2)*l/4 + rhoFluid*pi*(d^2)*l/4)/2;
%--------
% Analytical deflexion
delta = (1/2 * rhoFluid*cD0*d*va*va*l*l*l*l)/(8*E*I);
% structure damping
cu = 0;
% van der pol equation parameters
if ILVIVBool
    epsilony = 0.04; epsilonx = 0.02; St = 0.2; 
    B = 2 * pi * St * va / d;
    cq = epsilony * B; kq = B^2; 
    Ay = 12; Cy = Ay/d; %Cy = Ay/d;
    Ax = 96; Cx = Ax/d; %Cx = Ax/d 
else
    epsilony = 0.3; St = 0.2; 
    B = 2 * pi * St * va / d;
    cq = epsilony * B; kq = B^2; 
    Ay = 12; C2 = Ay/d; 
end
%--------
%
% Initial Conditions
%
q0 = 2; dq0 = 0; 
p0 = 2; dp0 = 0;
X0 = 0; dX0 = 0; ddX0 = 0; % 7.9763e-06
Y0 = 0; dY0 = 0; ddY0 = 0;
Z0 = 0; dZ0 = 0; ddZ0 = 0;
%--------
%
% ODE option
%options = odeset('RelTol',1e-6);
options = odeset('RelTol',1e-15);
%
