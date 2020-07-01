% Title: lum=PHluminosity
%
% Arguments: Bins (an array, preferrably row, can be column)
%            Alpha (The power the frequency goes)
% Returns: lum (number of new photons (a row vector))
% 
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 06/06/2020

function lum=PHluminosity(bins,alpha)
  k=max(size(bins));
  lum=zeros(k,1);
  lum=1e31.*(bins(1)./bins).^(alpha);
endfunction
