% Title: L=TEeHe
%
% Arguments: nh (Number density of Hydrogen (Neutral or Ionised) in m^-3 in each cell)
%            nhe (Number density of Helium (Neutral or Ionised) in m^-3 in each cell)
%            xh1 (Neutral hydrogen fraction in each cell) 
%            xhe1 (Neutral helium fraction in each cell)
%            xhe3 (Ionised helium fraction in each cell)
%            Temp (Temperature in each cell)
%
% Returns: L (Cooling rate from collisional excitation)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 26/08/2020

function L=TEeHe(nh,nhe,xh1,xhe1,xhe3,Temp)
  ne=(1.-xh1)*nh+(1.-xhe1+xhe3)*nhe;
  L=5.54e-30*ne.*nhe.*(1-xhe1-xhe3).*Temp.^(-0.397).*exp(-473638./Temp);
end