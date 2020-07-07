% Title: aT=REalphaII
%
% Arguments: Temp (Temperature, can be an array)
% Returns: aT (Recombination coefficient in m^3 s^-1)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 17/06/2020
%   Made Temperature array compatible 02/07/2020

function aT=REalphaII(Temp)
  aT=2.065e-17*sqrt(1./Temp).*(6.414.-0.5*log(Temp)+8.68e-3*(Temp.^(1/3)));
endfunction
