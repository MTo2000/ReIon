% Title: L=TEeH
%
% Arguments: nh (Number density of Hydrogen in m^-3) 
%            xh1 (Neutral Fraction (an array)) 
%            Temp (Temperature, can be an array)  
% Returns: L (Cooling rate from collisional excitation)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 30/06/2020

function L=TEeH(nh,nhe,xh1,xhe1,xhe3,Temp)
  ne=(1.-xh1)*nh+(1.-xhe1+xhe3)*nhe;
  L=7.3e-32*ne.*nh.*xh1.*e.^(-118400./Temp);
endfunction