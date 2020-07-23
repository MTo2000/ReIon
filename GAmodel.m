% Title: A=GAmodel
%
% Arguments: N (number matrix of photons before the effects of the gas)
%            x (array recording the position of each photon packets)
%            nh (Number density of hydrogens (neutral of ionised) in cm^-3) 
%            sigma (Cross section of hydrogens in cm^2) 
%            Nc (Integer number of cells)
%            L (Length of the box in cm)
%            xh1 (Neutral Fraction) 
%            bins (An array)
%            vT (Threshold frequency)
% Returns: A (number matrix of photons after the effects of the gas)
%          minxh (minimum neutral fractions to prevent numerical error growing exponentially)
%          total1 (total number of photons in each cell before the gas acts on it)
%          total2 (total number of photons in each cell after the gas acts on it)
%          totalG1 (total amount of photon energy that could be absorbed)
%          totalG2 (total amount of photon energy that exits) 
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 09/06/2020
%   Make tau an array 12/06/20
%   Added minimum calculations 25/06.20
%   Added total number calculations 28/06/20

function [A,minxh,totalH1,totalHe1,totalHe2,totalGH1,totalGHe1,totalGHe2]=GAmodel(N,x,nh,nhe,sigmaH1,sigmaHe1,sigmaHe2,Nc,L,xh1,xhe1,xhe3,bins,vH1,vHe1,vHe2)
  format short e
  
  xhe2=1.-xhe1-xhe3;
  
  tau1=L*nh*xh1'*sigmaH1/Nc;
  tau2=L*nhe*xhe1'*sigmaHe1/Nc;
  tau3=L*nhe*xhe2'*sigmaHe2/Nc;
  Tau=tau1+tau2+tau3;
  
  A=N;
  AH1=N.*(bins-vH1)*6.62607004e-34;
  AHe1=N.*(bins-vHe1)*6.62607004e-34;
  AHe2=N.*(bins-vHe2)*6.62607004e-34;
  
  mask1=find(AH1<0);
  mask2=find(AHe1<0);
  mask3=find(AHe2<0);
  AH1(mask1)=0;
  AHe1(mask2)=0;
  AHe2(mask3)=0;
  
  totalH1=zeros(1,Nc);
  totalHe1=zeros(1,Nc);
  totalHe2=zeros(1,Nc);
  totalGH1=zeros(1,Nc);
  totalGHe1=zeros(1,Nc);
  totalGHe2=zeros(1,Nc);
  
  pH1=e.^-tau1;
  pHe1=e.^-tau2;
  pHe2=e.^-tau3;
  qH1=1.-pH1;
  qHe1=1.-pHe1;
  qHe2=1.-pHe2;
  D=qH1.*pHe1.*pHe2+qHe1.*pHe2.*pH1+qHe2.*pH1.*pHe1;
  P1=qH1.*pHe1.*pHe2.*(1.-e.^(-Tau))./D;
  P2=qHe1.*pHe2.*pH1.*(1.-e.^(-Tau))./D;
  P3=qHe2.*pH1.*pHe1.*(1.-e.^(-Tau))./D;
  
  for k=1:Nc
    front=find(x>=(k-1)*L/Nc);
    back=find(x<k*L/Nc);
    aim=intersect(front,back);
    totalH1(k)=sum(sum(A(aim,:).*P1(k,:)));
    totalHe1(k)=sum(sum(A(aim,:).*P2(k,:)));
    totalHe2(k)=sum(sum(A(aim,:).*P3(k,:)));
    A(aim,:)=A(aim,:).*e.^(-Tau)(k,:);
    totalGH1(k)=sum(sum(AH1(aim,:).*P1(k,:)));
    totalGHe1(k)=sum(sum(AHe1(aim,:).*P2(k,:)));
    totalGHe2(k)=sum(sum(AHe2(aim,:).*P3(k,:)));
  endfor
  minxh=1e-10*Nc./(L*[nh*max(sigmaH1),nhe*max(sigmaHe1),nhe*max(sigmaHe2)]);
endfunction