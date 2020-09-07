% Title: L=TEphoton
%
% Arguments: nh (Number density of Hydrogen (Neutral or Ionised) in m^-3 in each cell)
%            nhe (Number density of Helium (Neutral or Ionised) in m^-3 in
%            each cell)
%            xh1 (Neutral hydrogen fraction in each cell)
%            xhe1 (Neutral helium fraction in each cell)
%            xhe3 (Ionised helium fraction in each cell) 
%            Temp (Temperature in each cell)  
%
% Returns: L (Cooling rate from photon cooling)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 30/06/2020

function L=TEphoton(nh,nhe,xh1,xhe1,xhe3,Temp)
  ne=(1.-xh1)*nh+(1.-xhe1+xhe3)*nhe;
  L=ne.*(REbetaHII(Temp).*(1-xh1)*nh+(REbetaHeII(Temp).*(1-xhe1-xhe3)+REbetaHeIII(Temp).*xhe3)*nhe);
end
