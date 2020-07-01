% Title: A=GAmodel
%
% Arguments: N (number matrix of photons before the effects of the gas)
%            x (array recording the position of each photon packets)
%            nh (Number density of hydrogens (neutral of ionised) in cm^-3) 
%            sigma (Cross section of hydrogens in cm^2) 
%            Nc (Integer number of cells)
%            L (Length of the box in cm)
%            xh1 (Neutral Fraction) 
% Returns: A (number matrix of photons after the effects of the gas)
%          minxh (minimum neutral fractions to prevent numerical error growing exponentially)
%          total1 (total number of photons in each cell before the gas acts on it)
%          total2 (total number of photons in each cell after the gas acts on it)
% 
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 09/06/2020
%   Make tau an array 12/06/20
%   Added minimum calculations 25/06.20
%   Added total number calculations 28/06/20

function [A,minxh,total1,total2]=GAmodel(N,x,nh,sigma,Nc,L,xh1)
  format short e
  tau=L*nh*xh1'*sigma/Nc;
  mask=find(tau<1e-10);
  A=N;
  total1=zeros(1,Nc);
  total2=zeros(1,Nc);
  for k=1:Nc
    front=find(x>=(k-1)*L/Nc);
    back=find(x<k*L/Nc);
    aim=intersect(front,back);
    total1(k)=sum(sum(A(aim,:)));
    A(aim,:)=A(aim,:).*e.^(-tau)(k,:);
    total2(k)=sum(sum(A(aim,:)));
  endfor
  minxh=1e-10*Nc/(L*nh*sigma(1));
endfunction