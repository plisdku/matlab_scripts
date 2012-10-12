function hh2he = matrixHH2HE(omega, k, x, varargin)
% matrixHH2HE Return the 2x2 matrix converting H_forward and H_backward
%   amplitudes and phases (measured from zero) to H and E at a point.
% 
% Usage: hh2he = matrixHH2HE(omega, k, x, epsr)
%   omega   angular frequency
%   k       longitudinal wavenumber, possibly complex
%   x       longitudinal position
%   epsr     (optional) relative permittivity
%
% The fields are assumed to be in SI units.  If k has a positive imaginary
% part, the forward fields will be attenuated.
%
% The inverse function is matrixHE2HH.

epsr = 1;

if nargin > 3
    epsr = varargin{1};
end

hh2he = [ exp(1i*k*x), exp(-1i*k*x); ...
          -k/omega/epsr*exp(1i*k*x), k/omega/epsr*exp(-1i*k*x) ];
