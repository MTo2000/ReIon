% Title: aT=REalphaHeII
%
% Arguments: Temp (Temperature, can be an array)
% Returns: aT (Recombination coefficient in m^3 s^-1)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 09/07/2020

function aT=REalphaHeII(Temp)
  aT=3.294e-17*(sqrt(Temp./15.54).*(1.+sqrt(Temp./15.54)).^(0.309).*(1.+sqrt(Temp./3.676e+7)).^(1.691)).^(-1).+1.9e-9.*(1.+0.3*e.^(-9.4e+4./Temp)).*e.^(-4.7e+5./Temp).*Temp.^(-1.5);
endfunction
