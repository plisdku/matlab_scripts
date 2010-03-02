% createDrudeMaterial   Set up a DrudeMetalModel material for Trogdor
%   Usage: createDrudeMaterial(name, epsinf, omegap, tauc)
%
%   Example:
%
%   goldMaterial = createDrudeMaterial('Gold', 12.99, 4e15, 9e-15);
%   b = createBlock([-5 -5 -5 5 5 5], 'Gold');
%   a = createAssembly(b);
%   g = createGrid(a, goldMaterial);
%
%   See also: createMaterial, createStaticDielectric, createPEC, createStaticLossyDielectric
%
%   version 4.5
%   July 29, 2008
function drudeMat = createDrudeMaterial(name, epsinf, omegap, tauc)

drudeMat = createMaterial(name, 'DrudeMetalModel', {'epsinf', num2str(epsinf)}, ...
    {'omegap', num2str(omegap)}, {'tauc', num2str(tauc)});