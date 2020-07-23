% Title: bT=REbetaHeII
%
% Arguments: Temp (Temperature, can be an array)
% Returns: bT (Recombination coefficient in J m^3 s^-1)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 09/07/2020

function bT=REbetaHeII(Temp)
  bT=1.55e-39*Temp.^(0.3647)+1.24e-26*(1.+0.3*e.^(-9.4e+4./Temp)).*e.^(-4.7e+5./Temp).*Temp.^(-1.5);
endfunction