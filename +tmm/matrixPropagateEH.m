function m = matrixPropagateEH(omega, k, distance, varargin)
% matrixPropagateEH Propagation matrix for E and H over a fixed distance
% 
% Usage: m = matrixPropagateEH(omega, k, distance, mur)
%   omega       angular frequency [1/s]
%   k           longitudinal wavenumber [1/m]
%   distance    longitudinal position [m]
%   mur         (optional) relative permeability [unitless]
%
% The fields are assumed to be in SI units.
%
% This matrix is equal to
%   matrixEE2EH(omega, k, distance, mur) * matrixEH2EE(omega, k, 0, mur)

mu0 = 4e-7*pi;
mur = 1;

if nargin > 3
    mur = varargin{1};
end

m = [ cos(k*distance), 1i*omega*mur*mu0/k*sin(k*distance); ...
      1i*k/omega/mur/mu0*sin(k*distance), cos(k*distance) ];

