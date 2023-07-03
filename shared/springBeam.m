function Udot = springBeam(t, y,m,ku,cu,kq,cq,l,d,rhof,cL0,cD0,cDi,Cx,Cy,va, qnp1k)
global VIVBool;
Udot = zeros( size(y) ) ;
% No q...
if length(y) == 2

  % y = [u, v]
  Uz = y(1) ;
  dUz = y(2) ;
  q =  2 ;
  Udot(1) = dUz ;
  Udot(2) = -ku/ m * Uz + 3*l/(8*m) * ( 1/2*rhof * cL0 * d * q/2 * ( va^2 + dUz^2 ) ) - cu*dUz ;
elseif length(y) == 4
  % y = [u, v, q, qdot]
  Uz = y(1) ;
  dUz = y(2) ;
  q = y(3) ;
  qdot = y(4) ;
  
  Udot(1) = dUz ;
  Udot(2) = -ku/ m * Uz + 3*l/(8*m) * ( 1/2*rhof * cL0 * d * q/2 * ( va^2 + dUz^2 ) ) - cu*dUz ;
  Udot(3) = qdot ;
  Udot(4) = - cq * (q^2-1)*qdot - kq * q + Cy * Udot(2) ;
%
elseif length(y) == 5
  % Lift and drag Without WOM
  % y = [Uz; dUz; X0; dX0; q0]
  Uz = y(1) ; % Lift
  dUz = y(2) ;
  Ux = y(3) ; % Drag
  dUx = y(4) ;
  qnp1k = y(5) ;
  %
  bl = 1/2*rhof * d * (cL0*qnp1k/2) * ( (-va-dUx)^2 + dUz^2 ) ;
  bd = 1/2*rhof * d * cD0 * ( (-va-dUx)^2 + dUz^2 ) ;
  %
  Udot(1) = dUz ;
  Udot(2) = -ku/ m * Uz + 3*l*bl/(8*m) - cu*dUz ;
  Udot(3) = dUx ;
  Udot(4) = -ku/ m * Ux - 3*l*bd/(8*m) - cu*dUx ; % Only one sign different
  Udot(5) = 0  ;

elseif length(y) == 6
    %y
  % Lift and drag
  % y = [Ux; dUx; X0; dX0; q0; dq0]
  Uz = y(1) ; % Lift
  dUz = y(2) ;
  Ux = y(3) ; % Drag
  dUx = y(4) ;
  q = y(5) ;
  qdot = y(6) ;
  %
  %q = 2; qdot =0;% For no VIV
  
%   bl = 1/2*rhof * d * (cL0*q/2) * ( (-va-dUx/2)^2 + (dUz/2)^2 ) ;
%   bd = 1/2*rhof * d * cD0 * ((-va-dUx/2)^2 + (dUz/2)^2 ) ;
%+ va if wind along ex
  bl = 1/2*rhof * d * (cL0*q/2) * ( (va-dUx/2)^2 + (dUz/2)^2 ) ;
  bd = 1/2*rhof * d * cD0 * ((va-dUx/2)^2 + (dUz/2)^2 ) ;
  
  Udot(1) = dUz ;
  Udot(2) = -ku/ m * Uz + 3*l*bl/(8*m) - cu*dUz ;
  Udot(3) = dUx ;
  Udot(4) = -ku/ m * Ux + 3*l*bd/(8*m) - cu*dUx ; % Only one sign different
  Udot(5) = qdot ;
  Udot(6) = - cq * (q^2-1)*qdot - kq * q + Cy * Udot(2) ;
  
elseif length(y) == 8
    % Lift and drag and ILVIV
  % y = [Ux; dUx; X0; dX0; q0; dq0; p0; dp0]
  Uz = y(1) ; % Lift
  dUz = y(2) ;
  Ux = y(3) ; % Drag
  dUx = y(4) ;
  q = y(5) ;
  qdot = y(6) ;
  p = y(7) ;
  pdot = y(8) ;
  %
  %q = 2; qdot =0;% For no VIV
  
%   bl = 1/2*rhof * d * (cL0*q/2) * ( (-va-dUx/2)^2 + (dUz/2)^2 ) ;
%   bd = 1/2*rhof * d * cD0 * ((-va-dUx/2)^2 + (dUz/2)^2 ) ;
%+ va if wind along ex
  bl = 1/2*rhof * d * (cL0*q/2) * ( (va-dUx/2)^2 + (dUz/2)^2 ) ;
  bd = 1/2*rhof * d * (cD0 + cDi*p/2)* ((va-dUx/2)^2 + (dUz/2)^2 ) ;
%   bl = 1/2*rhof * d * (cL0*q/2) * ( (va-dUx)^2 + (dUz)^2 ) ;
%   bd = 1/2*rhof * d * (cD0 + cDi*p/2)* ((va-dUx)^2 + (dUz)^2 ) ;
  
  Udot(1) = dUz ;
  Udot(2) = -ku/ m * Uz + 3*l*bl/(8*m) - cu*dUz ;
  Udot(3) = dUx ;
  Udot(4) = -ku/ m * Ux + 3*l*bd/(8*m) - cu*dUx ; % Only one sign different
  Udot(5) = qdot ;
  Udot(6) = - cq * (q^2-1)*qdot - kq * q + Cy * Udot(2)/2 ; % Average of nodal accelerations
  Udot(7) = pdot ;
  Udot(8) = - 2* cq/2 * (p^2-1)*pdot - 4 * kq * p + Cx * Udot(4)/2 ;
end