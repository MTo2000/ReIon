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

function A=PHsource(N,lum,bins,alpha,dt,rs,Ac)
  A=[dt.*lum(bins,alpha).*Ac./(4*pi*6.62607004e-27*rs^2);N];
endfunction
