% Title: A=GAmodel
%
% Arguments: N (number matrix of photons before the effects of the gas)
%            x (array recording the position of each photon packets)
%            ne (Number density of electrons in cm^-3) 
%            sigma (Cross section of electrons in cm^2) 
%            Nc (Integer number of cells)
%            L (Length of the box in cm)
%            Bins (An array) ('')
% Returns: A (number matrix of photons after the effects of the gas)
% 
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 09/06/2020
%   Make tau an array 12/06/20

function A=GAmode1(N,x,nh,sigma,Nc,L,xh1)
  tau=L*nh*xh1'*sigma/Nc
  b=size(x)(1);
  if b>Nc
    N(end-Nc+1:end,:)=N(end-Nc+1:end,:).*e.^(-tau);
    A=N;
  elseif
    A=N.*e.^(-flip(tau))(1:b,:);
  endif
endfunction