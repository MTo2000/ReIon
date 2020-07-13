% Title: L=TEphoton
%
% Arguments: nh (Number density of Hydrogen in m^-3) 
%            xh1 (Neutral Fraction (an array))
%            Recombination (Temperature dependent function) 
%            Temp (Temperature, float)  
% Returns: L (Cooling rate from photon cooling)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 30/06/2020

function L=TEphoton(nh,nhe,xh1,xhe1,xhe3,Temp)
  ne=(1.-xh1)*nh+(1.-xhe1+xhe3)*nhe;
  L=ne.*(REbetaHII(Temp).*(1.-xh1)*nh+(REbetaHeII(Temp).*(1-xhe1-xhe3)+REbetaHeIII(Temp).*xhe3)*nhe);
endfunction
