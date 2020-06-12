% Title: A=PHsource
%
% Arguments: N (number matrix before emission)
%            lum (Luminosity function)
%            Bins (An array)
%            Alpha (The power the frequency goes)
%            dt (Timestep)
% Returns: A (number matrix of photons after emission)
% 
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 06/06/2020

function A=PHsource(N,lum,bins,alpha,dt)
  A=[N;dt.*lum(bins,alpha)];
endfunction
