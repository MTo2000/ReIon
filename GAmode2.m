% Title: A=GAmode2
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

function A=GAmode2(N,x,nh,sigma,Nc,L,xh1)
  tau=L*nh*xh1'*sigma/Nc
  A=N;
  ibox=find(x<L);
  if length(ibox)>length(xh1)(1)
    error("Something is wrong with the code");
  endif
  A(ibox,:)=A(ibox,:).*e.^(-tau)(ibox,:);
endfunction