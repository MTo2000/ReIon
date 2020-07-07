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

function L=TEphoton(nh,xh1,Recombination,Temp)
  ne=(1.-xh1)*nh;
  L=Recombination(Temp).*ne.*ne;
endfunction
