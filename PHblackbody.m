% Title: lum=PHblackbody
%
% Arguments: Bins (an array, preferrably row, can be column)
%            T (Temperature, an array)
%            c (Speed of Light)
%            rad (Radius of source)
% Returns: lum (Luminosity)
% 
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 10/08/2020

function lum=PHblackbody(bins,T,c,rad)
  k=max(size(bins));
  lum=zeros(k,1);
  lum=(8/c^2)*6.62607004e-34*((pi*rad)^2)*(bins.^3)./(exp(6.62607004e-34*bins/(1.38064852e-23*T))-1);
end