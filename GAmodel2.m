% Title: A=GAmodel2
%
% Arguments: N (number matrix of photons before the effects of the gas)
%            x (array recording the position of each photon packets)
%            nh (Number density of hydrogens (neutral or ionised) in m^-3) 
%            nhe (Number density of heliums (neutral or ionised) in m^-3)
%            sigmaH1 (Array of Cross section of hydrogen I in m^2 for each frequency)
%            sigmaHe1 (Array of Cross section of helium I in m^2 for each frequency)
%            sigmaHe2 (Array of Cross section of helium II in m^2 for each frequency)
%            Nc (Integer number of cells)
%            L (Length of the box in m)
%            xh1 (Array of Neutral hydrogen fractions in each cell)
%            xhe1 (Array of Neutral helium fractions in each cell)
%            xhe3 (Ionised helium fraction in each cell)
%            bins (An array of frequency used)
%            vH1 (Threshold frequency for hydrogen I)
%            vHe1 (Threshold frequency for helium I)
%            vHe2 (Threshold frequency for helium II)
%            w (Weighting for each frequency, an array)
%
% Returns: A (number matrix of photons after the effects of the gas)
%          totalH1 (total number of hydrogen I that got excited in each cell)
%          totalHe1 (total number of helium I that got excited in each cell)
%          totalHe2 (total number of helium II that got excited in each cell)
%          totalGH1 (total energy released from hydrogen I excitation in each cell)
%          totalGHe1 (total energy released from helium I excitation in each cell)
%          totalGHe2 (total energy released from helium II excitation in
%          each cell) 
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 17/08/2020

function [A,totalH1,totalHe1,totalHe2,totalGH1,totalGHe1,totalGHe2]=GAmodel2(N,nh,nhe,sigmaH1,sigmaHe1,sigmaHe2,Nc,L,xh1,xhe1,xhe3,bins,vH1,vHe1,vHe2,w)
  format short e
  
  xhe2=1.-xhe1-xhe3;
  
  tau1=L*nh*xh1'*sigmaH1/Nc;
  tau2=L*nhe*xhe1'*sigmaHe1/Nc;
  tau3=L*nhe*xhe2'*sigmaHe2/Nc;
  Tau=tau1+tau2+tau3;
  
  %display(Tau(1:10,:))
  
  [a,b]=size(N);
  
  if a<Nc
    A=[N;zeros(Nc-a,b)];
  else
    A=N;
  end
  AH1=A.*(bins-vH1)*6.62607004e-34;
  AHe1=A.*(bins-vHe1)*6.62607004e-34;
  AHe2=A.*(bins-vHe2)*6.62607004e-34;
  
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
  
  pH1=exp(-tau1);
  pHe1=exp(-tau2);
  pHe2=exp(-tau3);
  qH1=1.-pH1;
  qHe1=1.-pHe1;
  qHe2=1.-pHe2;
  D=qH1.*pHe1.*pHe2+qHe1.*pHe2.*pH1+qHe2.*pH1.*pHe1;
  P1=qH1.*pHe1.*pHe2.*(1.-exp(-Tau))./D;
  P2=qHe1.*pHe2.*pH1.*(1.-exp(-Tau))./D;
  P3=qHe2.*pH1.*pHe1.*(1.-exp(-Tau))./D;

  totalH1=sum(A(1:Nc,:).*P1.*w,2)';
  totalHe1=sum(A(1:Nc,:).*P2.*w,2)';
  totalHe2=sum(A(1:Nc,:).*P3.*w,2)';
  A(1:Nc,:)=A(1:Nc,:).*exp(-Tau);
  totalGH1=sum(AH1(1:Nc,:).*P1.*w,2)';
  totalGHe1=sum(AHe1(1:Nc,:).*P2.*w,2)';
  totalGHe2=sum(AHe2(1:Nc,:).*P3.*w,2)';
  
  A=A(1:a,:);
end