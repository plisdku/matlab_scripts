function trogdor_end(varargin)
%trogdor_end End Trogdor simulation description and write parameter file
%   trogdor_end should be the last function called in describing a Trogdor 6
%   simulation.  It will write the parameter file params.xml and all other
%   data files necessary to run the given simulation.
%
%   trogdor_end('XML', paramFileName) permits customization of the
%   parameter file name.
%
%   trogdor_end('Directory', 'dat', 'OutputDirectory', 'outputs') will put all
%   the auxiliary files into the directory "dat" and all the field outputs
%   into the directory "outputs".

import com.mathworks.xml.XMLUtils.*;

X.XML = 'params.xml';
X.Directory = 'sim';
X.OutputDirectory = 'output';
X.Parameters = [];
X = parseargs(X, varargin{:});

if ~isstr(X.XML); error('Invalid filename'); end
if ~isstr(X.Directory); error('Invalid directory name'); end
if ~isstr(X.OutputDirectory); error('Invalid output directory name'); end

if ~isempty(X.Directory) && ~exist([pwd filesep X.Directory], 'dir')
    try mkdir([pwd filesep X.Directory])
    catch exception
        error('Could not create helper directory!');
    end
end

%{
if ~isempty(X.OutputDirectory) && ~exist([pwd filesep X.OutputDirectory], 'dir')
    try mkdir([pwd filesep X.OutputDirectory])
    catch exception
        error('Could not create output directory!');
    end
end
%}

% Last things:
%   store the extent of source grids in links

sim = t6.simulation();
sim.Directory = X.Directory;
sim.OutputDirectory = X.OutputDirectory;

% Make meshes if we're asked to
if ~isempty(sim.Grid.NodeGroup)
    sim.Grid.Meshes = sim.Grid.NodeGroup.meshes(X.Parameters);
end

doc = t6.xml.generateXML(sim);
xmlwrite([sim().directoryString, X.XML], doc);

delete(sim);

%t6.TrogdorSimulation.clear();
