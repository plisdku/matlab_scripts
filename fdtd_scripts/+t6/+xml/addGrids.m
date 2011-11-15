function addGrids(documentNode, simHandle)
doc = documentNode;
root = doc.getDocumentElement;
sim = simHandle;

for gg = 1:length(sim.Grids)
    grid = sim.Grids{gg};
    
    yeeCells = grid.YeeCells;
    nonPMLYeeCells(1:3) = yeeCells(1:3) + grid.PML(1:3);
    nonPMLYeeCells(4:6) = yeeCells(4:6) - grid.PML(4:6);
    %yeeCells(1:3) = nonPMLYeeCells(1:3) - grid.PML(1:3);
    %yeeCells(4:6) = nonPMLYeeCells(4:6) + grid.PML(4:6);
    
    if any(round(yeeCells) ~= yeeCells)
        error('Grid does not have integer dimensions.');
    end
    
    gridXML = doc.createElement('Grid');
    gridXML.setAttribute('name', grid.Name);
    gridXML.setAttribute('yeeCells', sprintf('%i ', yeeCells));
    gridXML.setAttribute('nonPML', sprintf('%i ', nonPMLYeeCells));
    gridXML.setAttribute('origin', sprintf('%i ', grid.Origin));
    %gridXML.setAttribute('nonPML',...
    %    sprintf('%i ', nonPMLYeeCells + [originTrogdor, originTrogdor]));
    %gridXML.setAttribute('origin', sprintf('%i ', originTrogdor));
    
    t6.xml.addAssembly(grid.Assembly, gridXML, doc, grid.Origin);
    t6.xml.addOutputs(grid, gridXML, doc, grid.Origin);
    t6.xml.addTFSFSources(grid, gridXML, doc, grid.Origin);
    t6.xml.addCustomTFSFSources(grid, gridXML, doc, grid.Origin);
    t6.xml.addHardSources(grid, gridXML, doc, grid.Origin);
    t6.xml.addCurrentSources(grid, gridXML, doc, grid.Origin);
    t6.xml.addPMLParams(grid, gridXML, doc);
    
    root.appendChild(gridXML);
end