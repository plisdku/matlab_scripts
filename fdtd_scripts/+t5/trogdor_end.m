function trogdor_end(varargin)
%trogdor_end End Trogdor simulation description and write parameter file
%   trogdor_end should be the last function called in describing a Trogdor 5
%   simulation.  It will write the parameter file params.xml and all other
%   data files necessary to run the given simulation.
%
%   trogdor_end('XML', paramFileName) permits customization of the
%   parameter file name.
%
%   trogdor_end('Directory', 'dat') will put all the auxiliary files into
%   the directory dat.

import com.mathworks.xml.XMLUtils.*;

X.XML = 'params.xml';
X.Directory = '';
X = parseargs(X, varargin{:});

if ~isstr(X.XML); error('Invalid filename'); end
if ~isstr(X.Directory); error('Invalid directory name'); end

if ~exist(X.Directory, 'dir')
    try mkdir(X.Directory)
    catch exception
        error('Could not create helper directory!');
    end
end

% Last things:
%   store the extent of source grids in links

sim = t5.TrogdorSimulation.instance;
sim.Directory = X.Directory;

for gg = 1:length(sim.Grids)
    grid = sim.Grids{gg};
    
    for ll = 1:length(grid.Links)
        sourceGridIndex = t5.indexOf(grid.Links{ll}.sourceGrid, sim.Grids);
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

doc = t5.xml.generateXML(sim);

directory = t5.TrogdorSimulation.instance().directoryString;

xmlwrite([directory, X.XML], doc);


