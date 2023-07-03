function windVel = windVelLinear(x,t)
%global finalTime; global vwindMax; 
global velocities;
  windx = velocities(t); % 
%   windx = t * vwindMax / finalTime
  windVel = [windx 0 0]'; 
end