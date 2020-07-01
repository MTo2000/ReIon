% Title: L=TCeH
%
% Arguments: nh (Number density of Hydrogen in cm^-3) 
%            xh1 (Neutral Fraction (an array)) 
%            Temp (Temperature, float)  
% Returns: L (Cooling rate from photon cooling)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 30/06/2020

function L=TCeH(nh,xh1,Temp)
  ne=(1.-xh1)*nh;
  L=7.3e-19*ne.*nh.*xh1.*e^(-118400/Temp);
endfunction