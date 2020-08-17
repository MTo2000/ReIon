function [ x, w ] = gen_legendre_compute ( x1, x2, n)

%*****************************************************************************80
%
%% GEN_LEGENDRE_COMPUTE computes a generalized Gauss-Legendre rule.
%
%  Discussion:
%
%    The integral to be approximated has the form:
%
%      Integral ( x1 < x < x2 ) f(x) dx
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    05 October 2012
%
%  Author:
%
%    Avery Meiksin
%
%  Reference:
%
%    W.H. Press, B. P. Flannery, S.A. Teukolsky, W.T. Vetterling
%    Numerical Recipes
%    Cambridge University Press, 1986,
%    ISBN: 0521308119,
%    LC: QA297.519.4.
%
%    Based on an algorithm due to G. Rybicki.
%
%  Parameters:
%
%    Input, integer x1, the lower x bound.
%
%    Input, integer x2, the upper x bound.
%
%    Input, integer n, the order of the rule.
%
%    Output, real x(n), w(n), the abscissas and weights
%    for the requested generalized Gauss-Legendre rule.
%
%
  EPS = 3.e-14;
  m = (n + 1)/ 2;
  xm = 0.5*(x2+x1);
  xl = 0.5*(x2-x1);
  
  z = cos(pi*(1 - .25)/ (n + .5));
  p1 = 1.;
  p2 = 0.;
    for j = 1:n
      p3 = p2;
      p2 = p1;
      p1 = ((2.*j - 1.)*z*p2 - (j - 1.)*p3)/ j;
    end
    pp = n*(z*p1 - p2)/ (z*z - 1.);
    z1 = z+p1/pp;
    
  for i = 1:m
    z = cos(pi*(i - .25)/ (n + .5));
    while abs(z-z1)>EPS
      p1 = 1.;
      p2 = 0.;
      for j = 1:n
        p3 = p2;
        p2 = p1;
        p1 = ((2.*j - 1.)*z*p2 - (j - 1.)*p3)/ j;
      end
      pp = n*(z*p1 - p2)/ (z*z - 1.);
      z1 = z;
      z = z1 - p1/pp;
    end
    x(i) = xm - xl*z;
    x(n + 1 - i) = xm + xl*z;
    w(i) = 2.*xl/ ((1. - z*z)*pp*pp);
    w(n + 1 - i) = w(i);
  end
  return
end
