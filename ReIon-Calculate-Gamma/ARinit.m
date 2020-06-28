% Initialise package of functions for computing atomic rates.
%
%DEFINE GLOBALS
%
%Physical constants
%
global m_p;
m_p = 1.6726216e-24;    %proton mass in g (CODATA 2006)
global kB;
kB = 1.3806504e-16;     %Boltzmann constant (erg/ K)
global hP;
hP = 6.6260690e-27;     %Planck constant (erg/ Hz)
global A10;
A10 = 2.85e-15;         %H n=1 hfs line transition rate in sec^{-1}
global lam10;
lam10 = 21.10611407;    %H n=1 hfs line wavelength in cm
global Y_He;
Y_He = 0.235;           %cosmic helium abundance
global T21cm;
T21cm = 0.06816864;     %21cm transition temperature (K)
global erg_per_eV;
erg_per_eV = 1.602176487e-12;  %erg/eV (CODATA 2006)
global Hz_per_eV;
Hz_per_eV = 2.417989454e14;    %Hz/eV (CODATA 2006)
%
%Set flags
%
global VERBOSE;
VERBOSE = false;
