% Title: A=PHbox
%
% Arguments: N_0 (number matrix at t=0)
%            dt (Timestep)
%            Gas (Function that describes the gas' effect on the photons)
%            L (Length of box)
%            lum (Luminosity function)
%            Bins (An array)
%            Alpha (The power the frequency goes)
% Returns: N (number matrix of photons after the first batch of photons reached the end of the box)
% 
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 04/06/2020

function [N,x]=PHbox(N_0,dt,gas,L,lum,bins,alpha)
  T=0;
  N=N_0;
  c=299792458;
  x=zeros(1,0);
  while c*T<L
    N=gas(N,c*T);
    T=T+dt;
    N=PHsource(N,lum,bins,alpha,dt);    #Source adding in photons
    x=c*dt.+x;
    x=[x;0];
  endwhile
%  display(N)
%  display(x)addasdsadasdada
endfunction