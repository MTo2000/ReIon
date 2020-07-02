% Title: T=TEtemp
%
% Arguments: S (Entropy Parameter)
%            xh1 (Neutral Fraction (an array))
%            rho (Mass density)  
% Returns: T (Temperature)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 02/07/2020

function T=TEtemp(S,rho,xh1)
  mb=1.6735575e-27./(2.-xh1);
  T=mb.*S*rho^(2/3)/1.38064852e-23;
endfunction
