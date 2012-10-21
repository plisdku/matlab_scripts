function addModalSource(varargin)
% addModalSource  Add electric and magnetic source currents to impress a field
%
% addModalSource('Bounds', [0 0 0 100 100 0], 'PhasorE', EE, 'PhasorH', HH, ...
%   'TimeEnvelope', @(t) exp(1i*freq*t), ...
%   'Direction', [0 0 1], 'X', xSrc, 'Y', ySrc, 'Z', zSrc)
% Drive the grid with an impressed field distribution (e.g. waveguide mode).
%
%   Usage: addModalSource(named parameters)
%
%   Named parameters:
%       Bounds      Region of simulation space in which to add electromagnetic
%                   current
%       PhasorE     An array of size [N M P 3] with complex amplitude of the E-
%                   field in a volume including Bounds.  Internally the fields
%                   necessary for the source currents will be interpolated to
%                   all necessary positions, so the number of samples N, M and P
%                   are at the discretion of the caller.
%       PhasorH     An array of size [N M P 3] with complex amplitude of the H-
%                   field.  These points should be colocated with E.
%       Direction   The "forward" direction of the mode
%       X           An N-element array of distances along the x-axis for each
%                   element of PhasorE or PhasorH
%       Y           An M-element array (see X)
%       Z           A P-element array (see X)
%       TimeFunction    A function of time returning a complex amplitude to
%                   multiply PhasorE and PhasorH by
%       Overwrite   "true" to rewrite the source J and M on every call to 
%                   trogdor_end.  "false" to save and re-use current sources.
%                   Setting Overwrite to false may save a lot of time setting
%                   up simulations that share source conditions!
%       Mode        "forward" or "adjoint".  If unspecified, the source will
%                   be added to both forward and adjoint simulations.
%
% Typically a modal source will be injected along a plane, or several planes,
% in the simulation space.  PhasorE and PhasorH should include samples
% surrounding the plane to guarantee successful field interpolation.  In
% particular, no dimension N, M or P of PhasorE or PhasorH should be less than
% two!  That is, we require that
%
%       X(1) < Bounds(1) <= Bounds(4) < X(end)
%       Y(1) < Bounds(2) <= Bounds(5) < Y(end)
%       Z(1) < Bounds(3) <= Bounds(6) < Z(end)
%

import t6.*
import modeInjection.*

X.Bounds = [];
X.Overwrite = true;
X.PhasorE = [];
X.PhasorH = [];
X.Direction = [];
X.X = [];
X.Y = [];
X.Z = [];
X.TimeFunction = [];
X.Mode = '';
X = parseargs(X, varargin{:});


[phasorJ phasorM] = modeInjection.surfaceEquivalentCurrents(...
    X.PhasorE, X.PhasorH, X.Direction, X.X, X.Y, X.Z);

srcJx = @(x,y,z) spatialModulation(phasorJ(x,y,z,1), X.TimeFunction);
srcJy = @(x,y,z) spatialModulation(phasorJ(x,y,z,2), X.TimeFunction);
srcJz = @(x,y,z) spatialModulation(phasorJ(x,y,z,3), X.TimeFunction);
srcMx = @(x,y,z) spatialModulation(phasorM(x,y,z,1), X.TimeFunction);
srcMy = @(x,y,z) spatialModulation(phasorM(x,y,z,2), X.TimeFunction);
srcMz = @(x,y,z) spatialModulation(phasorM(x,y,z,3), X.TimeFunction);

addCurrentSource('Field', 'jx jy jz mx my mz', 'Bounds', X.Bounds, ...
    'FieldFunctor', {srcJx srcJy srcJz srcMx srcMy srcMz}, ...
    'Overwrite', X.Overwrite, 'Mode', X.Mode);


function func = spatialModulation(modulationFunction, funcOfTime)
%spatialModulation  Create a function f(t) that returns a separable
% function of space and time such that f(t) = real(T(t)*X(x,y,z)).
%
% Usage:
%
% [xx yy zz] = ndgrid(xs, ys, zs);
% f = spatialModulation(sin(xx).*cos(yy), @(t) sin(t))
%
% This returns a function f that maps times t to arrays f(t), where
% size(f(t)) == size(xx) == size(yy) == size(zz).
%
% This function exists to speed up the process of writing current source
% data to files for Trogdor 6.  Often evaluating a function f(x,y,z,t) is
% more expensive than storing a modulating function M(x,y,z) and
% multiplying that data by a time-dependent signal T(t).
%

%func = @(t) real(modulationFunction funcOfTime(t));
func = @(t) bsxfun(@times, modulationFunction, ...
    reshape(funcOfTime(t), 1, 1, 1, []));