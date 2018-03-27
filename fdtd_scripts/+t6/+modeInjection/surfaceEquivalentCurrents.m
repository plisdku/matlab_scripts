function [phasorJ, phasorM] = surfaceEquivalentCurrents(srcE, srcH, inwardNormal, x0, y0, z0)

import t6.modeInjection.*

nCross = crossProductMatrix(inwardNormal);

% This means nCross*srcE where srcE is NxMx3 and nCross is 3x3.  The order
% of arguments to multTensor puts the tensor first, but it's really
% matrix*tensor!
srcM = t6.multTensor(srcE, nCross, 4);
srcJ = -t6.multTensor(srcH, nCross, 4);

% Neglecting wavevector still!  :-O
phasorM = @(x,y,z,xyz) gridInterp(x0,y0,z0,srcM(:,:,:,xyz),x,y,z);
phasorJ = @(x,y,z,xyz) gridInterp(x0,y0,z0,srcJ(:,:,:,xyz),x,y,z);
%phasorM = @(x,y,z,xyz) interpn(x0, y0, z0, srcM(:,:,:,xyz), x, y, z);
%phasorJ = @(x,y,z,xyz) interpn(x0, y0, z0, srcJ(:,:,:,xyz), x, y, z);

%{
% might be 10x slower this way
%phasorM = @(x,y,z,xyz) exp(1i*y*k_xy(2)) .* exp(1i*x*k_xy(1)) .* ...
%    reshape(spline(z_unitless, srcM(xyz,:), z(:)), size(z));
%phasorJ = @(x,y,z,xyz) exp(1i*y*k_xy(2)) .* exp(1i*x*k_xy(1)) .* ...
%    reshape(spline(z_unitless, srcJ(xyz,:), z(:)), size(z));

phasorM = @(x,y,z,xyz) exp(1i*y*k_xy(2)) .* exp(1i*x*k_xy(1)) .* ...
    reshape(interp1(z_unitless, srcM(xyz,:), z(:)), size(z));
phasorJ = @(x,y,z,xyz) exp(1i*y*k_xy(2)) .* exp(1i*x*k_xy(1)) .* ...
    reshape(interp1(z_unitless, srcJ(xyz,:), z(:)), size(z));
%}

