% Title: T=TEtemp
%
% Arguments: S (Entropy Parameter)
%            rho (Mass Density in kg m^-3)
%            xh1 (Neutral Hydrogen Fraction (an array))
%            xhe1 (Neutral Helium Fraction (an array))
%            xhe3 (Ionised Helium Fraction (an array))
%            R (Ratio)
% Returns: T (Temperature)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 02/07/2020

function T=TEtemp(S,rho,xh1,xhe1,xhe3,R)
  mb=(1+4/R)*1.6735575e-27./(2.+(2/R).-xh1-xhe1/R+xhe3/R);
  T=mb.*S*rho^(2/3)/1.38064852e-23;
endfunction
