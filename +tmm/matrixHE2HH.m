function he2hh = matrixHE2HH(omega, k, x, varargin)
%matrixHE2HH Return the 2x2 matrix converting H and E at a point to
%   H_forward and H_backward amplitudes and phases (measured from zero)
% 
% Usage: he2hh = matrixHE2HH(omega, k, x, epsr)
%   omega   angular frequency
%   k       longitudinal wavenumber, possibly complex
%   x       longitudinal position
%   epsr     (optional) relative permittivity
%
% The fields are assumed to be in SI units.
%
% The inverse function is matrixHH2HE.

epsr = 1;

if nargin > 3
    epsr = varargin{1};
end

he2hh = 0.5 * ...
    [ exp(-1i*k*x), -omega*epsr/k*exp(-1i*k*x); ...
      exp(1i*k*x), omega*epsr/k*exp(1i*k*x) ];
