% Title: A=GAmodel
%
% Arguments: N (number matrix of photons before the effects of the gas)
%            x (array recording the position of each photon packets)
%            ne (Number density of electrons in cm^-3) 
%            sigma (Cross section of electrons in cm^2) 
%            Nc (Integer number of cells)
%            L (Length of the box in cm)
% Returns: A (number matrix of photons after the effects of the gas)
% 
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 09/06/2020

function A=GAmode1(N,x,nh,sigma,Nc,L)
  tau=L*nh*sigma/Nc;
  b=size(x)(1);
  for a=1:b
    if x(a)<L
      A(a,:)=e^(-tau).*N(a,:);
    else
      A(a,:)=N(a,:);
    endif
  endfor
endfunction