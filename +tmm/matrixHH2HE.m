function hh2he = matrixHH2HE(omega, k, x, varargin)
% matrixHH2HE Return the 2x2 matrix converting H_forward and H_backward
%   amplitudes and phases (measured from zero) to H and E at a point.
% 
% Usage: ee2eh = matrixHH2HE(omega, k, x, mur)
%   omega   angular frequency [1/s]
%   k       longitudinal wavenumber, possibly complex [1/m]
%   x       longitudinal position [m]
%   epsr     (optional) relative permittivity [unitless]
%
% The fields are assumed to be in SI units.
%
% The inverse function is matrixHE2HH.

eps0 = 8.854187817e-12;
epsr = 1;

if nargin > 3
    epsr = varargin{1};
end

hh2he = [ exp(1i*k*x), exp(-1i*k*x); ...
          -k/omega/epsr/eps0*exp(1i*k*x), k/omega/epsr/eps0*exp(-1i*k*x) ];
