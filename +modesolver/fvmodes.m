function [ex ey ez hx hy hz neff] = fvmodes(lambda, guess, nmodes, dx, dy, ...
    varargin)
% fvmodes   Compute the electromagnetic field and effective index for
% waveguides.
%
% USAGE:
%
% [Ex Ey Ez Hx Hy Hz neff] = fvmodes(lambda, guess, nmodes, dx, dy, ...
%   eps, boundary)
% [Ex Ey Ez Hx Hy Hz neff] = fvmodes(lambda, guess, nmodes, dx, dy, ...
%   epsxx, epsyy, epszz, boundary)
% [Ex Ey Ez Hx Hy Hz neff] = fvmodes(lambda, guess, nmodes, dx, dy, ...
%   epsxx, epsxy, epsyx, epsyy, epszz, boundary)
%
%
% INPUT:
% 
% lambda - optical wavelength
% guess - scalar shift to apply when calculating the eigenvalues.
%     This routine will return the eigenpairs which have an
%     effective index closest to this guess
% nmodes - the number of modes to calculate
% dx - horizontal grid spacing (vector or scalar)
% dy - vertical grid spacing (vector or scalar)
% eps - index mesh (isotropic materials)  OR:
% epsxx, epsxy, epsyx, epsyy, epszz - index mesh (anisotropic)
% boundary - 4 letter string specifying boundary conditions to be
% applied at the edges of the computation window.  
%   boundary(1) = North boundary condition
%   boundary(2) = South boundary condition
%   boundary(3) = East boundary condition
%   boundary(4) = West boundary condition
% The following boundary conditions are supported: 
%   'A' - Hx is antisymmetric, Hy is symmetric.
%   'S' - Hx is symmetric and, Hy is antisymmetric.
%   '0' - Hx and Hy are zero immediately outside of the
%         boundary. 
% 
% OUTPUT:
% 
% Ex, Ey, Ez - calculated electric field, computed at the center of each
%   element
% Hx, Hy, Hz - calculated magnetic field, computed at the edges or vertices
%   of each element
%
% Hz - calculated longitudinal magnetic field.  This output will 
%   have the same dimensions as Hx and Hy.
% Ex, Ey, Ez - calculated electric field.  These field components 
%   are computed at the center of each element instead of on the
%   edges or vertices.
%
% NOTES:
%
% 1) The magnetic field components (Hx, Hy, and Hz) are
% calculated at the edges of each cell, whereas the electric
% field components are computed at the center of each cell.
% Therefore if size(eps) = [n,m], then the magnetic fields
% will have a size of [n+1,m+1] while the computed electric
% fields will have a size of [n,m].
%
% AUTHORS:  Thomas E. Murphy (tem@umd.edu)
%
% The convenience function fvmodes.m was added by
%
%           Paul C. Hansen (pch@stanford.edu)
%

[hx hy neff] = modesolver.wgmodes(lambda, guess, nmodes, dx, dy, varargin{:});

[hz ex ey ez] = modesolver.tools.postprocess(lambda, neff, ...
    hx, hy, dx, dy, varargin{:});

