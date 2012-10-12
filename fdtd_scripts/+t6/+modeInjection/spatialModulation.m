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