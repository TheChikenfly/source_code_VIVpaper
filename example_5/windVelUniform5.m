function windVel = windVelUniform5(x,t)
  global vwindMax;
%   finalTime = 0.2 ;
%   constWindTime = finalTime / 10 ;
  %vwindMax = 0.1 ;
  % constant profile
  windx =  vwindMax  ;
  % ramp profile
  %windx = t * vwindMax / constWindTime * (t <= constWindTime) + vwindMax * (t > constWindTime) ;
  windVel = [windx 0 0]'; 
end