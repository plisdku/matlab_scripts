function trogdor_end(varargin)
%trogdor_end End Trogdor simulation description and write parameter file
%   trogdor_end should be the last function called in describing a Trogdor 5
%   simulation.  It will write the parameter file params.xml and all other
%   data files necessary to run the given simulation.
%
%   trogdor_end(paramFileName) permits customization of the parameter file name.   

import com.mathworks.xml.XMLUtils.*;

paramFileName = 'params.xml';
if length(varargin) > 0
    if ~isstr(varargin{1})
        error('Invalid filename.');
    end
    paramFileName = varargin{:};
end

% Last things:
%   store the extent of source grids in links

sim = t6.TrogdorSimulation.instance;

for gg = 1:length(sim.Grids)
    grid = sim.Grids{gg};
    
    for ll = 1:length(grid.Links)
        sourceGridIndex = t6.indexOf(grid.Links{ll}.sourceGrid, sim.Grids);
        if sourceGridIndex == -1
            error('Unknown source grid %s', grid.Links{ll}.sourceGrid);
        end
        % compare this calculation of the origin in the source grid with the
        % calculation of the origin in xml.addGrids.m.
        srcNonPMLYeeCells = sim.Grids{sourceGridIndex}.extent;
        srcPML = sim.Grids{sourceGridIndex}.PML;
        srcTrogdorOrigin = -(srcNonPMLYeeCells(1:3) - srcPML(1:3));
        
        grid.Links{ll}.sourceTrogdorOrigin = srcTrogdorOrigin;
    end
end

doc = t6.xml.generateXML(sim);
xmlwrite(paramFileName, doc);
