% createStaticDielectric  Initialize dielectric material for Trogdor simulation
%   createStaticDielectric(name, epsr, mur) makes a Material object to be
%   parented by one or more Grids.  This is a convenience wrapper function for
%   createMaterial.
%
%   Example:
%
%   air = createStaticDielectric('Air', 1.0, 1.0);
%
%   g = createGrid('Main Grid', ..., air);
%
%   See also: createMaterial, createStaticLossyDielectric, createDrudeMaterial, createPEC
%
%   version 4.5
%   July 29, 2008
function sdMat = createStaticDielectric(name, epsr, mur)

sdMat = createMaterial(name, 'StaticDielectricModel', {'epsr', num2str(epsr)},...
    {'mur', num2str(mur)});