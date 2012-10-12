function m = matrixPropagateEH(omega, k, distance, varargin)
% matrixPropagateEH Propagation matrix for E and H over a fixed distance
% 
% Usage: m = matrixPropagateEH(omega, k, distance, mur)
%   omega       angular frequency
%   k           longitudinal wavenumber
%   distance    longitudinal position
%   mur         (optional) relative permeability
%
% This matrix is equal to
%   matrixEE2EH(omega, k, distance, mur) * matrixEH2EE(omega, k, 0, mur)

mur = 1;

if nargin > 3
    mur = varargin{1};
end

m = [ cos(k*distance), 1i*omega*mur/k*sin(k*distance); ...
      1i*k/omega/mur*sin(k*distance), cos(k*distance) ];

