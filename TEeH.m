% Title: L=TEeH
%
% Arguments: nh (Number density of Hydrogen (Neutral or Ionised) in m^-3 in each cell)
%            nhe (Number density of Helium (Neutral or Ionised) in m^-3 in
%            each cell)
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
%   Created in 30/06/2020

function L=TEeH(nh,nhe,xh1,xhe1,xhe3,Temp)
  ne=(1.-xh1)*nh+(1.-xhe1+xhe3)*nhe;
  L=7.5e-32*ne.*nh.*xh1.*exp(-118348./Temp);
end