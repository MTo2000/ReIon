%
%***********************
%*  sigma = ARAnu_HeI(y)
%***********************
%***********************
% Function to compute photoionization cross-section for HeI
% (Units: cm**2) 
%
% ARGUMENTS
%  y            nu/ nu_T where nu_T = threshold frequency (may be array)
%
%
% RETURNS
%  sigma        Cross-section (cm^2) (may be array)
%
% COMPATIBILITY: Matlab, Octave
%
% REQUIREMENTS: 
%  Initialise with ARinit.m.
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  11 09 15 Creation date.
%           (adapted from ARAnu_HZ.m)
% From Verner et al. (1996) ApJ 465:487-498
%
function sigma = ARAnu_HeI(y);
pi2 = 2*pi;
sig0 = 949.2e-18;
Eth = 24.59;
E0 = 13.61;
ya = 1.469;
yw = 2.039;
y0 = 0.4434;
y1 = 2.136;
P = 3.188;
sigma = zeros(1,length(y));
x = zeros(1,length(y));
yy = zeros(1,length(y));
Fyy = zeros(1,length(y));
maskn = find(y < 0);
if(~isempty(maskn))
  printf('ARAnu_HeI: index = %d\n',maskn);
  printf('ARAnu_HeI: nu/ nu_T = %g\n',y(maskn));
end
maskp = find(y >= 1);
if(~isempty(maskp))
  x(maskp) = Eth*y(maskp)/ E0 - y0;
  yy = (x.*x + y1*y1).^0.5;
  Fyy = ((x - 1).^2 + yw*yw).*yy.^(0.5*P - 5.5);
  Fyy = Fyy./ (1 + (yy/ ya).^0.5).^P;
  sigma(maskp) = sig0*Fyy(maskp);
end
