function addCurrentSources(sim, grid, gridXML, doc, mode)
global TROG_XML_COUNT___;

directory = sim.directoryString;

for ss = 1:length(grid.CurrentSources)
    src = grid.CurrentSources{ss};
    
    if ~isempty(src.mode) && ~strcmpi(src.mode, mode)
        continue;
    end
    
    t6.xml.addCurrentSource(src, sim, grid, gridXML, doc);
end




