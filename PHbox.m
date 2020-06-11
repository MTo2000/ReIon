% Title: A=PHbox
%
% Arguments: N_0 (number matrix at t=0)
%            Gas (Function that describes the gas' effect on the photons)
%            L (Length of box in cm)
%            lum (Luminosity function)
%            Bins (An array)
%            Alpha (The power the frequency goes)
%            ne (Number density of electrons in cm^-3) (')
%            sigma (Cross section of electrons in cm^2) (')
%            Testtime (Time the simulation will undergo) (')
%            nh (Number density of Hydrogen in cm^-3) (')
%            Ac (Cross section area of the box/cylinder in cm^-3)
% Returns: N (number matrix of photons after the first batch of photons reached the end of the box)
%          x (Location of each photon packets)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 04/06/2020
%   New Function 09/06/2020

function [N,x,Gamma]=PHbox(N_0,gas,L,lum,bins,alpha,ne,nh,Ac,sigma,iterations)
  recommend=10*(ne*sigma*L);
  display(["The recommended number of cells is around "  num2str(recommend)])
  Nc=input("Please enter the number of cells desired: ");
  if floor(Nc)==Nc                  #Check integer
    ;
  else
    error("Number of cells is not an integer")
  endif 
  T=0;
  N=N_0;
  c=299792458;
  dt=L/(Nc * c);
  x=zeros(1,1);
  tau=L*ne*sigma/Nc;
  Gamma=zeros(0,Nc);
  Nh=nh*Ac*L/Nc;
  while T<iterations
    k=size(x)(1);
    Gamma=[Gamma;zeros(1,Nc)];          #Add row to record Gamma for this iteration
    total=sum(N,2)';                  #Record total number of photons entering each cells
    N=gas(N,x,ne,sigma,Nc,L);         #Gas acts on the photons
    T=T+1;
    x=c*dt.+x;
    for n=1:min(k,Nc)
      Gamma(end,n)=total(end-n+1)*(1-e^(-tau))./(dt*Nh);
    endfor
    x=[x;0];
    N=PHsource(N,lum,bins,alpha,dt);    #Source adding in photons
  endwhile
%  display(N)
%  display(x)
endfunction