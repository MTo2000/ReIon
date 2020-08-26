% Title: N,x=PHbox
%
% Arguments: Gas (Function that describes the gas' effect on the photons)
%            L (Length of box in cm)
%            nh (Number density of Hydrogen in m^-3)
%            ratio (Ratio of Hyydrogen atoms to Helium atoms) 
%            Ac (Cross section area of the box/cylinder in m^3)
%            duration (How long will the simulation take place)
%            xh1i (Initial Neutral Hydrogen Fraction)
%            xhe1i (Initial Neutral Helium Fraction)
%            xhe3i (Initial Fully Ionised Helium Fraction)
%            T0 (Initial Temperature) 
%            rs (Some initial distances from source) 
%            QN (Number of nodes for Quadrature)
%            rad (Radius of source)
% Returns: N (number matrix of photons after the first batch of photons reached the end of the box)
%          x (Location of each photon packets)
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
%   Helium included 07/07/20
%   Exclusive to Blackbody tweek 10/08/2020
%   Special case when M==1 14/08/2020
%   Change of graph outputs 18/08/2020

function [N,x]=PHbox(L,nh,ratio,Ac,duration,xh1i,xhe1i,xhe3i,T0,rs,QN,rad)
  
  format long e
  
  if xhe1i<0 || xhe3i<0
    error ("Helium fractions can't be negative")
  end
  if xhe1i+xhe3i>1
    error("Fraction of neutral and fully ionised Helium can't exceed 1")
  end
  
  rho=nh*1.6735575e-27*(1+4/ratio); %Mass Density in kg m^-3
  c=299792458;    %Speed of light
  vH1=3.282e+15;
  vHe1=5.933e+15;
  vHe2=1.313e+16;
  
  [bin1,w1]=gen_legendre_compute(vH1,vHe1,QN);
  [bin2,w2]=gen_legendre_compute(vHe1,vHe2,QN);
  [u,w3]=gen_legendre_compute(0,1/vHe2,QN);
  
  bin3=flip(1./u);
  w1=w1./bin1;
  w2=w2./bin2;
  w3=flip(w3).*bin3;
  
  bins=[bin1,bin2,bin3];
  w=[w1,w2,w3];
  
  BinH1=bins./vH1;
  BinHe1=bins./vHe1;
  BinHe2=bins./vHe2;
  sigmaH1=ARAnu_HZ(BinH1,1)*1e-4; %Turning it to m^2
  sigmaHe1=ARAnu_HeI(BinHe1)*1e-4; %Turning it to m^2
  sigmaHe2=ARAnu_HZ(BinHe2,2)*1e-4; %Turning it to m^2
  a=max([max(sigmaH1),max(sigmaHe1),max(sigmaHe2)]);
  recommend=10*(nh*xh1i*a*L);
  display(["The recommended number of cells is around "  num2str(recommend)])
  Nc=input("Please enter the number of cells desired: ");
  if floor(Nc)==Nc                  %Check integer
  else
    error("Number of cells is not an integer")
  end
  if xh1i<=1 && xh1i>0
  else
    error("Neutral fraction must be between 0 and 1")
  end
  dt=L/(Nc * c);
  modifier=dt/2*1e-10; 
  display(["The unmodified timestep is " num2str(dt)]);
  display(["The recommended modifier is " num2str(modifier)]);
  M=input("Please enter the modifier desired: ");
  dt=dt/M;
  it=ceil(duration/dt)
  T=1
  x=zeros(1,0);
  s=size(bins);
  N=zeros(0,s(2));
  nhe=nh/ratio;
  Nh=nh*Ac*L/Nc;                       %Total number of gas in each cell
  Nhe=Nh/ratio;
  xh1=xh1i*ones(1,Nc);                 %Initial fractions for all cells
  xhe1=xhe1i*ones(1,Nc);
  xhe3=xhe3i*ones(1,Nc);
  Temp=T0*ones(1,Nc);
  S=TEentropy(Temp,rho,xh1,xhe1,xhe3,ratio);
  
  %go=fopen('Gamma_list.txt','w+');
  %x1o=fopen('xh1_list.txt','w+');
  %x2o=fopen('xhe1_list.txt','w+');
  %x3o=fopen('xhe3_list.txt','w+');
  %to=fopen('Temp_list.txt','w+');
  %Go=fopen('G_list.txt','w+');
  %Lo=fopen('L_list.txt','w+');
  
  g_list=zeros(1,it);
  l_list=zeros(1,it);
  t_list=zeros(1,it);
  T_list=T0*ones(1,it);
  
  if M==1
    while T<it
      x=[rs;x];
    
      GammaH1=zeros(1,Nc);
      GammaHe1=zeros(1,Nc);
      GammaHe2=zeros(1,Nc);
      GH1=zeros(1,Nc);
      GHe1=zeros(1,Nc);
      GHe2=zeros(1,Nc);
    
      N=PHsource(N,PHblackbody(bins,1e+5,c,rad),dt,rs,Ac);    %Source adding in photons
      %N=N.*(rs./x).^2;
      [N,mini,totalH1,totalHe1,totalHe2,totalGH1,totalGHe1,totalGHe2]=GAmodel2(N,nh,nhe,sigmaH1,sigmaHe1,sigmaHe2,Nc,L,xh1,xhe1,xhe3,bins,vH1,vHe1,vHe2,w);         %Gas acts on the photons
    
      %neg1=find(xh1<mini(1));
      %xh1(neg1)=mini(1);
      %neg2=find(xhe1<mini(2));
      %xhe1(neg2)=mini(2);
      %xhe2=1.-xhe1.-xhe3;
      %neg3=find(xhe2<mini(3));
      %xhe2(neg3)=mini(3);
      %xhe3=max(1.-xhe1.-xhe2,1e-10);
      nonz=find(totalH1~=0);
       
      xhe2=1-xhe1-xhe3;
      
      GammaH1(nonz)=totalH1(nonz)./(dt*Nh*xh1(nonz));
      GammaHe1(nonz)=totalHe1(nonz)./(dt*Nhe*xhe1(nonz));
      non2=find(xhe2~=0);
      non3=intersect(nonz,non2);
      GammaHe2(non3)=totalHe2(non3)./(dt*Nhe*xhe2(non3));
      GH1(nonz)=nh*totalGH1(nonz)./(dt*Nh);
      GHe1(nonz)=nhe*totalGHe1(nonz)./(dt*Nhe);
      GHe2(nonz)=nhe*totalGHe2(nonz)./(dt*Nhe);  

    
      l=TEeH(nh,nhe,xh1,xhe1,xhe3,Temp)+TEphoton(nh,nhe,xh1,xhe1,xhe3,Temp);
      ne=(1.-xh1)*nh+(1.-xhe1+xhe3)*nhe;
      xh1=xh1+dt*(-GammaH1.*xh1+REalphaHII(Temp).*(1.-xh1).*ne);
      xhe1=xhe1+dt*(-GammaHe1.*xhe1+REalphaHeII(Temp).*xhe2.*ne);
      xhe3=xhe3+dt*(GammaHe2.*xhe2-REalphaHeIII(Temp).*xhe3.*ne);
      G=GH1+GHe1+GHe2;
      S=S+dt*2/3*rho^(-5/3)*(G-l);
      Temp=TEtemp(S,rho,xh1,xhe1,xhe3,ratio);
      %neg1=find(xh1<mini(1));
      %xh1(neg1)=mini(1);
      %neg2=find(xhe1<mini(2));
      %xhe1(neg2)=mini(2);
      %xhe2=1.-xhe1.-xhe3;
      %neg3=find(xhe2<mini(3));
      %xhe2(neg3)=mini(3);
      %xhe3=1.-xhe1.-xhe2;
    
      x=c*dt+x;
    
      g_list(T+1)=G(1);
      l_list(T+1)=l(1);
      T_list(T+1)=Temp(1);
      t_list(T+1)=T*dt;
      
      T=T+1
      
      %fprintf(go,"Time = %2e     ",T*dt);
      %fprintf(go,"%2e  ",GammaH1);
      %fprintf(go,"\n");
      %fprintf(x1o,"Time = %2e     ",T*dt);
      %fprintf(x1o,"%2e  ",xh1);
      %fprintf(x1o,"\n");
      %fprintf(x2o,"Time = %2e     ",T*dt);
      %fprintf(x2o,"%2e  ",xhe1);
      %fprintf(x2o,"\n");
      %fprintf(x3o,"Time = %2e     ",T*dt);
      %fprintf(x3o,"%2e  ",xhe3);
      %fprintf(x3o,"\n");
      %fprintf(to,"Time = %2e     ",T*dt);
      %fprintf(to,"%2e  ",Temp);
      %fprintf(to,"\n");
      %fprintf(Go,"Time = %2e     ",T*dt);
      %fprintf(Go,"%2e  ",G);
      %fprintf(Go,"\n");
      %fprintf(Lo,"Time = %2e     ",T*dt);
      %fprintf(Lo,"%2e  ",l);
      %fprintf(Lo,"\n");
  
    end
  else
    while T<it
      x=[rs;x];
    
      GammaH1=zeros(1,Nc);
      GammaHe1=zeros(1,Nc);
      GammaHe2=zeros(1,Nc);
      GH1=zeros(1,Nc);
      GHe1=zeros(1,Nc);
      GHe2=zeros(1,Nc);
    
      N=PHsource(N,PHblackbody(bins,1e+5,c,rad),dt,rs,Ac);    %Source adding in photons
      N=N.*(rs./x).^2;
      [N,mini,totalH1,totalHe1,totalHe2,totalGH1,totalGHe1,totalGHe2]=GAmodel(N,x,nh,nhe,sigmaH1,sigmaHe1,sigmaHe2,Nc,L,xh1,xhe1,xhe3,bins,vH1,vHe1,vHe2,w,rs);         %Gas acts on the photons
    
      %neg1=find(xh1<mini(1));
      %xh1(neg1)=mini(1);
      %neg2=find(xhe1<mini(2));
      %xhe1(neg2)=mini(2);
      %xhe2=1.-xhe1.-xhe3;
      %neg3=find(xhe2<mini(3));
      %xhe2(neg3)=mini(3);
      %xhe3=max(1.-xhe1.-xhe2,1e-10);
    
      filled=floor(T/M);
      ext=rem(T,M);
      nonz=find(totalH1~=0);
    
      for z=nonz
        if z<=filled
          GammaH1(z)=totalH1(z)./(dt*Nh*xh1(z)*M);
          GammaHe1(z)=totalHe1(z)./(dt*Nhe*xhe1(z)*M);
          if 1-xhe1(z)-xhe3(z)~=0
            GammaHe2(z)=totalHe2(z)./(dt*Nhe*(1-xhe1(z)-xhe3(z))*M);
          else
            GammaHe2(z)=0;
          end
          GHe1(z)=nhe*totalGHe1(z)./(dt*Nhe*M);
          GHe2(z)=nhe*totalGHe2(z)./(dt*Nhe*M);
        else
          GammaH1(z)=totalH1(z)./(dt*Nh*xh1(z)*ext);
          GammaHe1(z)=totalHe1(z)./(dt*Nhe*xhe1(z)*ext);
          if 1-xhe1(z)-xhe3(z)~=0
            GammaHe2(z)=totalHe2(z)./(dt*Nhe*(1-xhe1(z)-xhe3(z))*ext);
          else
            GammaHe2(z)=0;
          end
          GH1(z)=nh*totalGH1(z)./(dt*Nh*ext);
          GHe1(z)=nhe*totalGHe1(z)./(dt*Nhe*ext);
          GHe2(z)=nhe*totalGHe2(z)./(dt*Nhe*ext);
        end
      end
    
      l=TEeH(nh,nhe,xh1,xhe1,xhe3,Temp)+TEphoton(nh,nhe,xh1,xhe1,xhe3,Temp);
      ne=(1.-xh1)*nh+(1.-xhe1+xhe3)*nhe;
      xh1=xh1+dt*(-GammaH1.*xh1+REalphaHII(Temp).*(1.-xh1).*ne);
      xhe1=xhe1+dt*(-GammaHe1.*xhe1+REalphaHeII(Temp).*(1.-xhe1-xhe3).*ne);
      xhe3=xhe3+dt*(GammaHe2.*(1-xhe1-xhe3)-REalphaHeIII(Temp).*xhe3.*ne);
      G=GH1+GHe1+GHe2;
      S=S+dt*2/3*rho^(-5/3)*(G-l);
      Temp=TEtemp(S,rho,xh1,xhe1,xhe3,ratio);

      %neg1=find(xh1<mini(1));
      %xh1(neg1)=mini(1);
      %neg2=find(xhe1<mini(2));
      %xhe1(neg2)=mini(2);
      %xhe2=1.-xhe1.-xhe3;
      %neg3=find(xhe2<mini(3));
      %xhe2(neg3)=mini(3);
      %xhe3=1.-xhe1.-xhe2;
    
      x=c*dt+x;
    
      g_list(T+1)=G(1);
      l_list(T+1)=l(1);
      t_list(T+1)=T*dt;
      
      T=T+1
      
      %fprintf(go,"Time = %2e     ",T*dt);
      %fprintf(go,"%2e  ",GammaH1);
      %fprintf(go,"\n");
      %fprintf(x1o,"Time = %2e     ",T*dt);
      %fprintf(x1o,"%2e  ",xh1);
      %fprintf(x1o,"\n");
      %fprintf(x2o,"Time = %2e     ",T*dt);
      %fprintf(x2o,"%2e  ",xhe1);
      %fprintf(x2o,"\n");
      %fprintf(x3o,"Time = %2e     ",T*dt);
      %fprintf(x3o,"%2e  ",xhe3);
      %fprintf(x3o,"\n");
      %fprintf(to,"Time = %2e     ",T*dt);
      %fprintf(to,"%2e  ",Temp);
      %fprintf(to,"\n");
      %fprintf(Go,"Time = %2e     ",T*dt);
      %fprintf(Go,"%2e  ",G);
      %fprintf(Go,"\n");
      %fprintf(Lo,"Time = %2e     ",T*dt);
      %fprintf(Lo,"%2e  ",l);
      %fprintf(Lo,"\n");
  
    end
  end
  
  %fclose(go);
  %fclose(x1o);
  %fclose(x2o);
  %fclose(x3o);
  %fclose(to);
  %fclose(Go);
  %fclose(Lo);

  ratio0=l_list./g_list;
  ratio1=l./G;
  
  xp=c*dt*[1:Nc]+rs;
  
  limt=t_list(end);
  limx=xp(end);
  
  subplot(2,4,1)
  plot(t_list,ratio0);
  axis([0,limt,0,1.2*max(ratio0)])
  xlabel("Time in seconds")
  title("Heating Rate vs Cooling Rate ratio")
  
  subplot(2,4,2)
  plot(t_list,g_list);
  hold on
  plot(t_list,l_list);
  hold off
  axis([0,limt,0,1.2*max(g_list)])
  xlabel("Time in seconds")
  ylabel("Energy change in Joules per second")
  title("Heating Rate vs Cooling Rate, global")
  
  subplot(2,4,3)
  plot(t_list,T_list);
  axis([0,limt,0,1.2*max(T_list)])
  xlabel("Time in seconds")
  ylabel("Temperature in K")
  title("Temperature vs Time graph (Nearest cell)")
  
  subplot(2,4,4)
  plot(xp,Temp);
  axis([rs,limx,0,1.2*max(Temp)])
  xlabel("Distance in meters")
  ylabel("Temperature in K")
  title("Temperature vs Distance graph")
  
  subplot(2,4,5)
  plot(xp,ratio1);
  axis([rs,limx,0,1.2*max(ratio1)])
  xlabel("Distance in meters")
  ylabel("Heating and Cooling ratio")
  title("Ratio vs Distance graph")
  
  subplot(2,4,6)
  plot(xp,xh1);
  axis([rs,limx,0,1.2])
  xlabel("Distance in meters")
  ylabel("Neutral Hydrogen Fraction")
  title("Neutral Hydrogen Fraction")
  
  subplot(2,4,7)
  plot(xp,xhe1);
  axis([rs,limx,0,1.2])
  xlabel("Distance in meters")
  ylabel("Neutral Helium Fraction")
  title("Neutral Helium Fraction")
  
  subplot(2,4,8)
  plot(xp,xhe3);
  axis([rs,limx,0,1.2])
  xlabel("Distance in meters")
  ylabel("Ionised Helium Fraction")
  title("Ionised Helium Fraction")
  
  T_list=log(T_list);
  T2_list=zeros(1,length(T_list)-1);
  h=length(T_list)-1;
  
  for f=1:h
    T2_list(f)=T_list(f+1)-T_list(f);
    T2_list(f)=T2_list(f)/dt;
  end
  
  T2_list=[1,T2_list];
  
  %figure(2)
  %plot(t_list,1./T2_list)
  
end