function documentNode = generateXML(simulationHandle)
sim = simulationHandle; % shorthand

% The automatically-generated source files are numbered sequentially.
% This counts them.
global TROG_XML_COUNT___;
TROG_XML_COUNT___ = struct(...
    'softTime', 0,...
    'softMask', 0,...
    'softSpaceTime', 0,...
    'hardTime', 0,...
    'hardMask', 0,...
    'hardSpaceTime', 0,...
    'currentTime',0,...
    'currentSpaceTime',0,...
    'tfsfTime',0,...
    'heightMap',0);

doc = com.mathworks.xml.XMLUtils.createDocument('Simulation');

root = doc.getDocumentElement;
root.setAttribute('version', '6.0');
root.setAttribute('dx', num2str(sim.Dxyz(1), 10));
root.setAttribute('dy', num2str(sim.Dxyz(2), 10));
root.setAttribute('dz', num2str(sim.Dxyz(3), 10));
root.setAttribute('dt', num2str(sim.Dt, 10));
root.setAttribute('numT', num2str(sim.NumT));

t6.xml.addMaterials(doc, sim);
t6.xml.addGrids(doc, sim);


documentNode = doc;

