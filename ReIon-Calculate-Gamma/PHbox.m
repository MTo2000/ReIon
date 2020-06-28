% Title: A=PHbox
%
% Arguments: Gas (Function that describes the gas' effect on the photons)
%            L (Length of box in cm)
%            lum (Luminosity function)
%            Bins (An array)
%            Alpha (The power the frequency goes)
%            ne (Number density of electrons in cm^-3) (')
%            sigma (Cross section of electrons in cm^2) (')
%            Testtime (Time the simulation will undergo) (')
%            nh (Number density of Hydrogen in cm^-3) (')
%            Ac (Cross section area of the box/cylinder in cm^-3)(')
%            vT (Threshold Frequency)('')
%            xh1i (Initial neutral Fraction)(''')
%            Recombination (Temperature dependent function) (''')
%            Temp (Temperature, float) (''')
%            rs (Some initial distances from source) ('''')
% Returns: N (number matrix of photons after the first batch of photons reached the end of the box)
%          x (Location of each photon packets)
%          xh (Neutral Hydrogen Fraction, an array) ('''')
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 04/06/2020
%   New Function 09/06/2020
%   Make Tau an array 12/06/2020
%   Write Gamma in a separate file 12/06/20
%   Neutral Fraction Calculation 16/06/20

function [N,x,xh1]=PHbox(gas,L,lum,bins,alpha,nh,Ac,iterations,vT,xh1i,recombination,Temp,rs)
  format short e
  T=0;
  c=299792458;
  Z=input("Please enter the atomic number of the gas: ");
  Bins=bins./vT;
  sigma=ARAnu_HZ(Bins,Z);
  recommend=10*(nh*xh1i*max(sigma)*L);
  display(["The recommended number of cells is around "  num2str(recommend)])
  Nc=input("Please enter the number of cells desired: ");
  if floor(Nc)==Nc                  #Check integer
    ;
  elseif
    error("Number of cells is not an integer")
  endif
  if xh1i<1 && xh1i>0
    ;
  elseif
    error("Neutral fraction must be between 0 and 1")
  endif
  dt=L/(Nc * c);
  x=zeros(1,0);
  N=zeros(0,size(bins)(2));
  Nh=nh*Ac*L/Nc                        #Total number of Hydrogens in each cell
  xh1=xh1i*ones(1,Nc);                 #Initial neutral fraction for all cells
  if isfile('Gamma_list.txt')
     delete 'Gamma_list.txt';
  else
     ;
  endif
  go=fopen('Gamma_list.txt','w');
  while T<iterations
    x=[0;x];
    total1=zeros(1,Nc);
    total2=zeros(1,Nc);
    Gamma=zeros(1,Nc);
    N=PHsource(N,lum,bins,alpha,dt,rs,Ac);    #Source adding in photons
    k=find(x<=(Nc-0.1)*L/Nc);
    total1(k)=sum(N,2)(k);
    [N,mini]=gas(N,x,nh,sigma,Nc,L,xh1);         #Gas acts on the photons
    total2(k)=sum(N,2)(k);
    T=T+1;
    x=c*dt.+x;
    neg=find(xh1<mini);
    xh1(neg)=mini;
    Gamma=(total1-total2)./(dt*Nh*xh1);
    fprintf(go,"Time = %f     ",T*dt);
    fprintf(go,"%2e  ",Gamma);
    fprintf(go,"\n");
    #xh1=NFhydrogen1(xh1,Gamma,recombination,Temp,nh);
    xh1=xh1+dt*(-Gamma.*xh1+(recombination(Temp)*nh).*(1.-xh1).^2);
    neg=find(xh1<mini);
    xh1(neg)=mini;
  endwhile
  fclose(go);
endfunction