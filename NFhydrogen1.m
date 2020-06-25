% Title: x=NFhydrogen1
%
% Arguments: xh1 (Neutral Fraction of neutral hydrogen for last iteration)
%            Gamma (Ionisation rate, an array)
%            Recombination (Temperature dependent function)
%            Temp (Temperature, float)
%            Nh (Total number of hudrogen atoms (regardless or ionisd or not))
% Returns: x1 (Neutral Fraction of neutral hydrogen for next iteration)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 17/06/2020

function x1=NFhydrogen1(xh1,gamma,recombination,Temp,nh)
  format short e
  x1=zeros(1,length(xh1));
  ne=(1.-xh1).*nh; 
  non=find(gamma~=0); #Gas reached the cell
  zer=find(gamma==0); #Gas not reached the cell
  T=recombination(Temp);
  x1(non)=1./(gamma(non)./(ne(non).*T).+1);
  x1(zer)=xh1(zer);
endfunction