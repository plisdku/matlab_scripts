% createPEC     Initialize PEC material for Trogdor simulation
%   createPEC(name) makes a Material object to be parented by one or more Grids.
%   This is a convenience wrapper function for createMaterial.
%
%   Example:
%
%   pecMat = createPEC('PEC');
%
%   g = createGrid('Main Grid', ..., pecMat);
%
%   See also: createMaterial, createDrudeMaterial, createStaticDielectric, createStaticLossyDielectric
%
%   version 4.5
%   July 29, 2008
function pecMat = createPEC(name)

pecMat = createMaterial(name, 'PECModel');