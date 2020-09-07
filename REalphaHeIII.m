% Title: aT=REalphaHeIII
%
% Arguments: Temp (Temperature, can be an array)
%
% Returns: aT (Recombination coefficient in m^3 s^-1)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 09/07/2020

function aT=REalphaHeIII(Temp)
  aT=8.26e-17*sqrt(1./Temp).*(7.107-0.5*log(Temp)+5.47e-3*(Temp.^(1/3)));
end
