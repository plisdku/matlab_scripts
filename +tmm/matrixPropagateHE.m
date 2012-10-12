function m = matrixPropagateHE(omega, k, distance, varargin)
% matrixPropagateHE Propagation matrix for H and E over a fixed distance
% 
% Usage: m = matrixPropagateHE(omega, k, distance, mur)
%   omega       angular frequency
%   k           longitudinal wavenumber
%   distance    longitudinal position
%   epsr         (optional) relative permeability [unitless]
%
% This matrix is equal to
%   matrixHH2HE(omega, k, distance, epsr) * matrixHE2HH(omega, k, 0, epsr)

epsr = 1;

if nargin > 3
    epsr = varargin{1};
end

m = [ cos(k*distance), -1i*omega*epsr/k*sin(k*distance); ...
      -1i*k/omega/epsr*sin(k*distance), cos(k*distance) ];

