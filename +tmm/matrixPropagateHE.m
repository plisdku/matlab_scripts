function m = matrixPropagateHE(omega, k, distance, varargin)
% matrixPropagateHE Propagation matrix for H and E over a fixed distance
% 
% Usage: m = matrixPropagateHE(omega, k, distance, mur)
%   omega       angular frequency [1/s]
%   k           longitudinal wavenumber [1/m]
%   distance    longitudinal position [m]
%   epsr         (optional) relative permeability [unitless]
%
% The fields are assumed to be in SI units.
%
% This matrix is equal to
%   matrixHH2HE(omega, k, distance, epsr) * matrixHE2HH(omega, k, 0, epsr)

eps0 = 8.854187817e-12;
epsr = 1;

if nargin > 3
    epsr = varargin{1};
end

m = [ cos(k*distance), -1i*omega*epsr*eps0/k*sin(k*distance); ...
      -1i*k/omega/epsr/eps0*sin(k*distance), cos(k*distance) ];

