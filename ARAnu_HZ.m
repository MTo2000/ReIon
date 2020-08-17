%
%************************
%*  sigma = ARAnu_HZ(y,Z)
%************************
%************************
% Function to compute photoionization cross-section for hydrogenic atoms
% of atomic number Z. (Units: cm**2) 
%
% ARGUMENTS
%  y            nu/ nu_T where nu_T = threshold frequency (may be array)
%  Z            Atomic number of ion (Z=1 for H+; Z = 2 for He++)
%
%
% RETURNS
%  sigma        Cross-section (cm^2) (may be array)
%
% COMPATIBILITY: Matlab(?), Octave
%
% REQUIREMENTS: 
%  Initialise with ARinit.m.
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  11 03 13 Creation date.
%           (adapted from anu_HZ.f)
%
function sigma = ARAnu_HZ(y,Z);
pi2 = 2*pi;
anu_HZ = zeros(1,length(y));
sigma = zeros(1,length(y));
anu_H = 6.304318e-18;    %HI nu_L photoelectric cross-section (cm^2) (CODAT06)
maskn = find(y < 0);
if(~isempty(maskn))
  printf('ARAnu_HZ: index = %d\n',maskn);
  printf('ARAnu_HZ: nu/ nu_T = %g\n',y(maskn));
end
maskp = find(y >= 1);
if(~isempty(maskp))
  eps(maskp) = (y(maskp) - 1.).^0.5;
  maskt = find(eps < 1.e-5);
  maskl = find(eps >= 1.e-5);
  maskhn = find(eps < 0.5);
  maskhp = find(eps >= 0.5);
  ii = intersect(maskhn,maskl);
  jj = intersect(maskhp,maskl);
  kk = intersect(maskp,maskt);
  if(~isempty(kk))
    sigma(kk) = anu_H/ Z/ Z;
  end
  if(~isempty(ii))
    anu_HZ(ii) = y(ii).^4./ exp(4.*(1. - atan(eps(ii))./ eps(ii)));
    sigma(ii) = (anu_H/ (Z*Z)).*(1./ anu_HZ(ii));
  end
  if(~isempty(jj))
    anu_HZ(jj) = 1. - 1./ exp(pi2./ eps(jj));
    anu_HZ(jj) = anu_HZ(jj).* y(jj).^4./ exp(4.*(1. - atan(eps(jj))./ eps(jj)));
    sigma(jj) = (anu_H/ (Z*Z)).*(1./ anu_HZ(jj));
  end
end
