function addOutputs(sim, grid, gridXML, doc, mode)

%directory = sim.outputDirectoryString;

for oo = 1:length(grid.Outputs)
    output = grid.Outputs{oo};
    
    if ~isempty(output.mode) && ~strcmpi(output.mode, mode)
        continue;
    end
    
    t6.xml.addOutput(output, sim, grid, gridXML, doc);
end
