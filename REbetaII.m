% Title: bT=REbetaII
%
% Arguments: Temp (Temperature, can be an array)
% Returns: bT (Recombination coefficient in J m^3 s^-1)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 30/06/2020
%   Made Temperature array compatible 02/07/2020

function bT=REbetaII(Temp)
  bT=2.851e-40*sqrt(Temp).*(5.914.-0.5*log(Temp)+0.01184.*(Temp.^(1/3)));
endfunction
