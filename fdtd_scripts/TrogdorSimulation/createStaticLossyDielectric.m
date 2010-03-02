% createStaticLossyDielectric  Initialize lossy dielectric material for Trogdor simulation
%   createStaticLossyDielectric(name, epsr, mur, sigma) makes a Material object
%   to be parented by one or more Grids.  This is a convenience wrapper function
%   for createMaterial.
%
%   Example:
%
%   air = createStaticLossyDielectric('Air', 1.0, 1.0, 0.0);
%
%   g = createGrid('Main Grid', ..., air);
%
%   If working with complex permittivity, sigma should be epsi/omega at a
%   frequency of interest.
%
%   See also: createMaterial, createStaticDielectric, createDrudeMaterial, createPEC
%
%   version 4.5
%   July 29, 2008
function sdMat = createStaticLossyDielectric(name, epsr, mur, sigma)

sdMat = createMaterial(name, 'StaticLossyDielectricModel', {'epsr', num2str(epsr)},...
    {'mur', num2str(mur)}, {'sigma', num2str(sigma)});