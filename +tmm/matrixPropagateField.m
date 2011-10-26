function m = matrixPropagateField(k, distance)
% matrixPropagateField Propagation matrix for forward and backward E or H
%   fields
% 
% Usage: m = matrixPropagate(k, distance)
%   k           longitudinal wavenumber [1/m]
%   distance    longitudinal position [m]
%

m = diag([exp(1i*k*distance), exp(-1i*k*distance)]);