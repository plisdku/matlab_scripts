function addModalSource(varargin)
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