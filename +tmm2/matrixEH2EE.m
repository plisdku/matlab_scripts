function eh2ee = matrixEH2EE(omega, k, x, varargin)
%matrixEH2EE Return the 2x2 matrix converting Ez and Hy at a point x to
%   E_forward and E_backward amplitudes and phases (measured from zero)
% 
% Usage: eh2ee = matrixEH2EE(omega, k, x, mur)
%   omega   angular frequency
%   k       longitudinal wavenumber, possibly complex
%   x       longitudinal position
%   mur     (optional) relative permeability
%
% The inverse function is matrixEE2EH.

mur = 1;

if nargin > 3
    mur = varargin{1};
end

eh2ee = 0.5 * ...
    [ exp(-1i*k*x), -omega*mur/k*exp(-1i*k*x); ...
      exp(1i*k*x), omega*mur/k*exp(1i*k*x) ];
