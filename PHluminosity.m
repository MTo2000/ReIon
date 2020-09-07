% Title: lum=PHluminosity
%
% Arguments: Bins (an array, preferrably row, can be column)
%            Alpha (The power the frequency goes)
%            vT (Threshold frequency)
%
% Returns: lum (Luminosity)
% 
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 06/06/2020

function lum=PHluminosity(bins,alpha,vT)
  k=max(size(bins));
  lum=zeros(k,1);
  lum=1e24.*(vT./bins).^(alpha);
end