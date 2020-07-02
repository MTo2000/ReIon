% Title: S=TEentropy
%
% Arguments: Temp (Temperature, an array)
%            xh1 (Neutral Fraction (an array))
%            rho (Mass density)  
% Returns: S (Entropy Parameter, an array)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 02/07/2020

function S=TEentropy(Temp,rho,xh1)
  mb=1.6735575e-27./(2.-xh1);
  S=1.38064852e-23*rho^(-2/3).*Temp./mb;
endfunction
