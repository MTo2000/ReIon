% Title: A=PHsource
%
% Arguments: N (number matrix before emission)
%            flux (energy flux of photons per second)
%            dt (Timestep)
%            rs (Starting distance)
%            Ac (Cross section of area of the box/cylinder in m^3)
%            vT (Threshold Frequency)
% Returns: A (number matrix of photons after emission)
% 
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 06/06/2020
%   Changed compatibility requirement 10/08/2020

function A=PHsource(N,flux,dt,rs,Ac,vT)
  A=[dt.*flux.*Ac./(4*pi*6.62607004e-34*rs^2);N];
endfunction