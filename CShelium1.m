% Title: sigma=CShelium1
%
% Arguments: y (nu/nu_T, where nu_T = threshold frequency) 
% Returns: sigma (Cross-section in m^2, can be an array)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 07/07/2020

function sigma=CShelium1(y)
  suf=find(y>=1);
  sigma=zeros(1,length(y));
  sigma(suf)=6.30e-22*(1.66*y(suf).^(-2.05)-0.66*y(suf).^(-3.05));
endfunction
