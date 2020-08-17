function lum=PHblackbody2(bins,T)
  k=max(size(bins));
  lum=zeros(k,1);
  lum=104.22*(bins.^2)./(e.^(6.62607004e-34*bins/(1.38064852e-23*T)).-1);
endfunction