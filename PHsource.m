% Title: A=PHsource
%
% Arguments: N (number matrix before emission)
%            flux (energy flux of photons per second)
%            dt (Timestep)
%            vT (Threshold Frequency)
%
% Returns: A (number matrix of photons after emission)
% 
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 06/06/2020
%   Changed compatibility requirement 10/08/2020

function A=PHsource(N,flux,v,dt)
  %A=[dt.*flux.*Ac./(4*pi*6.62607004e-34*rs^2);N];
  A=[dt.*flux./(6.62607004e-34*v);N];
end