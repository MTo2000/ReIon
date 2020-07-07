% Title: N,x,xh1,eq=PHbox
%
% Arguments: Gas (Function that describes the gas' effect on the photons)
%            L (Length of box in cm)
%            lum (Luminosity function)
%            bins (An array)
%            alpha (The power the frequency goes)
%            nh (Number density of Hydrogen in m^-3) 
%            Ac (Cross section area of the box/cylinder in m^3)
%            iterations (Number of iterations)
%            vT (Threshold Frequency)
%            xh1i (Initial Neutral Fraction)
%            T0 (Initial Temperature) 
%            rs (Some initial distances from source) 
% Returns: N (number matrix of photons after the first batch of photons reached the end of the box)
%          x (Location of each photon packets)
%          xh (Neutral Hydrogen Fraction, an array) 
%          eq (Equilibrium value for Neutral Hydrogen Fraction)
%
% Compatibility: Octave (+Matlab?)
% Author: To Kwok Hei Matthew
% History:
%   Created in 04/06/2020
%   New Function 09/06/2020
%   Make Tau an array 12/06/2020
%   Write Gamma in a separate file 12/06/20
%   Equilibirum Neutral Fraction Calculation 16/06/20
%   Time dependent Neutral Fraction Calculation 25/06/20
%   Write xh1 in a separate file 30/06/20
%   Temperature dependence added 01/07/20

function [N,x,xh1,eq]=PHbox(gas,L,lum,bins,alpha,nh,Ac,iterations,vT,xh1i,T0,rs)
  format short e
  rho=nh*1.6735575e-27; #Mass Density in kg m^-3
  T=1;
  c=299792458;
  Bins=bins./vT;
  sigmaH1=ARAnu_HZ(Bins,1).*1e-4; #Turning it to m^2
  recommend=10*(nh*xh1i*max(sigmaH1)*L);
  display(["The recommended number of cells is around "  num2str(recommend)])
  Nc=input("Please enter the number of cells desired: ");
  if floor(Nc)==Nc                  #Check integer
    ;
  elseif
    error("Number of cells is not an integer")
  endif
  if xh1i<=1 && xh1i>0
    ;
  elseif
    error("Neutral fraction must be between 0 and 1")
  endif
  dt=L/(Nc * c);
  modifier=dt/2*1e-10; 
  display(["The unmodified timestep is " num2str(dt)]);
  display(["The recommended modifier is " num2str(modifier)]);
  M=input("Please enter the modifier desired: ");
  dt=dt/M;
  eq=zeros(1,Nc);
  x=zeros(1,0);
  N=zeros(0,size(bins)(2));
  Nh=nh*Ac*L/Nc;                        #Total number of Hydrogens in each cell
  xh1=xh1i*ones(1,Nc);                 #Initial neutral fraction for all cells
  Temp=T0*ones(1,Nc);
  S=TEentropy(Temp,rho,xh1);
  if isfile('Gamma_list.txt')
     delete 'Gamma_list.txt';
  else
     ;
  endif
  if isfile('xh_list.txt')
     delete 'xh_list.txt';
  else
     ;
  endif
  if isfile('Temp_list.txt')
     delete 'Temp_list.txt';
  else
     ;
  endif
  if isfile('G_list.txt')
     delete 'G_list.txt';
  else
     ;
  endif
  if isfile('L_list.txt')
     delete 'L_list.txt';
  else
     ;
  endif
  go=fopen('Gamma_list.txt','w');
  xo=fopen('xh_list.txt','w');
  to=fopen('Temp_list.txt','w');
  Go=fopen('G_list.txt','w');
  Lo=fopen('L_list.txt','w');
  xh_list=[xh1i];
  g_list=[0];
  l_list=[0];
  t_list=[0];
  T_list=[T0];
  G_list=[0];
  R_list=[0];
  while T<iterations
    x=[0;x];
    Gamma=zeros(1,Nc);
    G=zeros(1,Nc);
    N=PHsource(N,lum,bins,alpha,dt,rs,Ac,vT);    #Source adding in photons
    [N,mini,total1,total2,totalG1,totalG2]=gas(N,x,nh,sigmaH1,Nc,L,xh1,bins,vT);         #Gas acts on the photons
    neg=find(xh1<mini);
    xh1(neg)=mini;
    filled=floor(T/M);
    ext=rem(T,M);
    nonz=find(total1~=0);
    for z=nonz
      if z<=filled
        Gamma(z)=(total1(z)-total2(z))./(dt*Nh*xh1(z)*M);
        G(z)=nh*(totalG1(z)-totalG2(z))./(dt*Nh*M);
      elseif
        Gamma(z)=(total1(z)-total2(z))./(dt*Nh*xh1(z)*ext);
        G(z)=nh*(totalG1(z)-totalG2(z))./(dt*Nh*ext);
      endif
    endfor
    l=TEeH(nh,xh1,Temp)+TEphoton(nh,xh1,@REbetaII,Temp);
    xh1=xh1+dt*(-Gamma.*xh1+(REalphaII(Temp)*nh).*(1.-xh1).^2);
    S=S+dt*2/3 *rho^(-5/3)*(G-l);
    Temp=TEtemp(S,rho,xh1);
    neg=find(xh1<mini);
    xh1(neg)=mini;
    T=T+1;
    x=c*dt.+x;
    xh_list=[xh_list,xh1(1)];
    g_list=[g_list,G(1)];
    l_list=[l_list,l(1)];
    T_list=[T_list,Temp(1)];
    t_list=[t_list,T*dt];
    G_list=[G_list,xh1(1).*Gamma(1)];
    R_list=[R_list,((REalphaII(Temp)*nh).*(1.-xh1(1)).^2)(1)];
    fprintf(go,"Time = %2e     ",T*dt);
    fprintf(go,"%2e  ",Gamma);
    fprintf(go,"\n");
    fprintf(xo,"Time = %2e     ",T*dt);
    fprintf(xo,"%2e  ",xh1);
    fprintf(xo,"\n");
    fprintf(to,"Time = %2e     ",T*dt);
    fprintf(to,"%2e  ",Temp);
    fprintf(to,"\n");
    fprintf(Go,"Time = %2e     ",T*dt);
    fprintf(Go,"%2e  ",G);
    fprintf(Go,"\n");
    fprintf(Lo,"Time = %2e     ",T*dt);
    fprintf(Lo,"%2e  ",l);
    fprintf(Lo,"\n");
  endwhile
  fclose(go);
  fclose(xo);
  fclose(to);
  fclose(Go);
  fclose(Lo);
  eq=NFhydrogen1(xh1,Gamma,@REalphaII,Temp,nh);
  xh_list=eq(1)./xh_list;
  
  limx=t_list(end);
  lim1=5*min(max(g_list),max(l_list));
  lim2=5*min(max(G_list),max(R_list));
  
  subplot(2,3,1)
  plot(t_list,xh_list);
  axis([0,limx,0,1.2*max(xh_list)])
  xlabel("Time in seconds")
  legend("eq(xh1)/xh1(t)")
  title("How quickly it is approaching equilibrium value")
  
  subplot(2,3,2)
  plot(t_list,g_list);
  hold on
  plot(t_list,l_list);
  hold off
  axis([0,limx,0,lim1])
  xlabel("Time in seconds")
  ylabel("Energy change in Joules per second")
  legend("Heating Rate","Cooling Rate")
  title("Heating Rate vs Cooling Rate, zoomed in")
  
  subplot(2,3,3)
  plot(t_list,g_list);
  hold on
  plot(t_list,l_list);
  hold off
  axis([0,limx])
  xlabel("Time in seconds")
  ylabel("Energy change in Joules per second")
  legend("Heating Rate","Cooling Rate")
  title("Heating Rate vs Cooling Rate, global")
  
  subplot(2,3,4)
  plot(t_list,T_list);
  axis([0,limx,0,1.2*max(T_list)])
  xlabel("Time in seconds")
  ylabel("Temperature in K")
  legend("Temperature(t)")
  title("Temperature")
  
  subplot(2,3,5)
  plot(t_list,G_list);
  hold on
  plot(t_list,R_list);
  hold off
  axis([0,limx,0,lim2])
  xlabel("Time in seconds")
  ylabel("Rate of Change")
  legend("Ionisation Rate","Recombination Rate")
  title("Ionisation Rate vs Recombination Rate, zoomed in")
  
  subplot(2,3,6)
  plot(t_list,G_list);
  hold on
  plot(t_list,R_list);
  hold off
  axis([0,limx])
  xlabel("Time in seconds")
  ylabel("Rate of Change")
  legend("Ionisation Rate","Recombination Rate")
  title("Ionisation Rate vs Recombination Rate, global")
  
endfunction