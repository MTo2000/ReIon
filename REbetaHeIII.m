% Title: bT=REbetaHeIII
%
% Arguments: Temp (Temperature, can be an array)
%
% Returns: bT (Recombination coefficient in J m^3 s^-1)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 09/07/2020

function bT=REbetaHeIII(Temp)
  bT=1.14e-39*sqrt(Temp).*(6.607-0.5*log(Temp)+7.459e-3.*(Temp.^(1/3)));
end