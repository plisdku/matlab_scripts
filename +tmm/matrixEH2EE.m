function eh2ee = matrixEH2EE(omega, k, x, varargin)
%matrixEH2EE Return the 2x2 matrix converting E and H at a point to
%   E_forward and E_backward amplitudes and phases (measured from zero)
% 
% Usage: eh2ee = matrixEH2EE(omega, k, x, mur)
%   omega   angular frequency [1/s]
%   k       longitudinal wavenumber, possibly complex [1/m]
%   x       longitudinal position [m]
%   mur     (optional) relative permeability [unitless]
%
% The fields are assumed to be in SI units.
%
% The inverse function is matrixEE2EH.

mu0 = 4e-7*pi;
mur = 1;

if nargin > 3
    mur = varargin{1};
end

eh2ee = 0.5 * ...
    [ exp(-1i*k*x), omega*mur*mu0/k*exp(-1i*k*x); ...
      exp(1i*k*x), -omega*mur*mu0/k*exp(1i*k*x) ];
