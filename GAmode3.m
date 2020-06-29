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

function [A,minxh,total1,total2]=GAmode3(N,x,nh,sigma,Nc,L,xh1)
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