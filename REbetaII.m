% Title: bT=REbetaII
%
% Arguments: Temp (Temperature, float)
% Returns: bT (Recombination coefficient in cm^3 s^-1)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 30/06/2020

function bT=REbetaII(Temp)
  bT=2.851e-34*sqrt(Temp)*(5.914-0.5*log(Temp)+0.01184*(Temp^(1/3)));
endfunction
