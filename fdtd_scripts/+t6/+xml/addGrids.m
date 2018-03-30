function addGrids(doc, sim, designParameters, mode)
root = doc.getDocumentElement;

%grid = sim.Grids{gg};

yeeCells = sim.Grid.YeeCells;
nonPMLYeeCells(1:3) = yeeCells(1:3) + sim.Grid.PML(1:3);
nonPMLYeeCells(4:6) = yeeCells(4:6) - sim.Grid.PML(4:6);

if any(round(yeeCells) ~= yeeCells)
    error('Grid does not have integer dimensions.');
end

gridXML = doc.createElement('Grid');
gridXML.setAttribute('name', sim.Grid.Name);
gridXML.setAttribute('yeeCells', sprintf('%i ', yeeCells));
gridXML.setAttribute('nonPML', sprintf('%i ', nonPMLYeeCells));
gridXML.setAttribute('origin', sprintf('%i ', sim.Grid.Origin));
%gridXML.setAttribute('nonPML',...
%    sprintf('%i ', nonPMLYeeCells + [originTrogdor, originTrogdor]));
%gridXML.setAttribute('origin', sprintf('%i ', originTrogdor));

t6.xml.addAssembly(sim, gridXML, doc, sim.Grid.Origin);
t6.xml.addOutputs(sim, gridXML, doc, mode);
t6.xml.addTFSFSources(sim, gridXML, doc, mode);
t6.xml.addCustomTFSFSources(sim, gridXML, doc, mode);
t6.xml.addHardSources(sim, gridXML, doc, mode);
t6.xml.addCurrentSources(sim, gridXML, doc, mode);
t6.xml.addPMLParams(sim, gridXML, doc);
t6.xml.addMeasurements(sim, gridXML, doc, mode);

root.appendChild(gridXML);
