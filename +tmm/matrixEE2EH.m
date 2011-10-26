function ee2eh = matrixEE2EH(omega, k, x, varargin)
% matrixEE2EH Return the 2x2 matrix converting E_forward and E_backward
%   amplitudes and phases (measured from zero) to E and H at a point.
% 
% Usage: ee2eh = matrixEE2EH(omega, k, x, mur)
%   omega   angular frequency [1/s]
%   k       longitudinal wavenumber, possibly complex [1/m]
%   x       longitudinal position [m]
%   mur     (optional) relative permeability [unitless]
%
% The fields are assumed to be in SI units.
%
% The inverse function is matrixEH2EE.

mu0 = 4e-7*pi;
mur = 1;

if nargin > 3
    mur = varargin{1};
end

ee2eh = [ exp(1i*k*x), exp(-1i*k*x); ...
          k/omega/mur/mu0*exp(1i*k*x), -k/omega/mur/mu0*exp(-1i*k*x) ];
