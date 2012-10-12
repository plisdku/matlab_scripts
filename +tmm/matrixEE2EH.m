function ee2eh = matrixEE2EH(omega, k, x, varargin)
% matrixEE2EH Return the 2x2 matrix converting E_forward and E_backward
%   amplitudes and phases (measured from zero) to E and H at a point.
% 
% Usage: ee2eh = matrixEE2EH(omega, k, x, mur)
%   omega   angular frequency
%   k       longitudinal wavenumber, possibly complex
%   x       longitudinal position
%   mur     (optional) relative permeability
%
% The fields are assumed to be in SI units.  If k has a positive imaginary
% part, the forward fields will be attenuated.
%
% The inverse function is matrixEH2EE.

%mu0 = 4e-7*pi;
mur = 1;

if nargin > 3
    mur = varargin{1};
end

ee2eh = [ exp(1i*k*x), exp(-1i*k*x); ...
          k/omega/mur*exp(1i*k*x), -k/omega/mur*exp(-1i*k*x) ];
